"""llm2geneset: using LLMs to generate gene sets."""

import asyncio
import re
from importlib import resources
from typing import List

import json_repair
import pandas as pd
import tiktoken
import tqdm.asyncio
from asynciolimiter import StrictLimiter
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def read_gmt(gmt_file: str):
    """
    Load a GMT file.

    See for details:
    https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

    Args:
       gmt_file: a gene set file in Broad GMT format
    Returns:
       (descr, genes): tuple of descriptions and list of genes
    """
    with open(gmt_file, "r") as file:
        lines = file.readlines()

    descr = []
    genes = []
    for line in lines:
        sp = line.strip().split("\t")
        sp2 = [re.sub(",.*", "", value) for value in sp[2:]]
        sp2 = [x.upper() for x in sp2 if x]
        if len(sp2) > 0:
            descr.append(sp[0])
            genes.append(sp2)

    return (descr, genes)


def get_embeddings(
    client, text_list: List[str], model="text-embedding-3-large", batchsz=2048
):
    """Get embeddings using OpenAI API, processing in batches of 2048.

    Args:
        client: synchronous OpenAI client
        text_list: lists of texts to embed
        model: embedding model
        batchsz: size of batches max is 2048
    Returns:
        List of embeddings.
    """

    def chunks(lst, n):
        """Yield successive n-sized chunks from lst (sliding window)."""
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    # Lowercase all text entries
    text_list_cleaned = []
    for text in text_list:
        text_list_cleaned.append(text.lower())

    # Process in batches of 2048
    all_embeddings = []
    for batch in chunks(text_list_cleaned, 2048):  # for each sliding window
        response = client.embeddings.create(model=model, input=batch).data
        # Extract embeddings for the current batch (sliding window) and append to
        # the all_embeddings list
        batch_embeddings = [resp.embedding for resp in response]
        all_embeddings.extend(batch_embeddings)
    # for the text_list, convert sliding windows of text into sliding
    # windows of embeddings
    return all_embeddings


def extract_last_code_block(markdown_text):
    """Extract last code block in output.

    Args:
       markdown_text: text with markdown blocks
    Returns:
       Returns last code block. Raises exception if no code block
       was found.
    """
    # Regular expression to find code blocks enclosed in triple backticks
    pattern = r"```(?:[\w]*\n)?([\s\S]*?)```"
    code_blocks = re.findall(pattern, markdown_text)
    if not code_blocks:
        raise ValueError("No code blocks found")
    # Return the last code block, if any
    return code_blocks[-1]


async def get_genes(
    aclient,
    descr,
    model="gpt-4o",
    prompt_type="basic",
    use_sysmsg=False,
    n_retry=3,
    use_tqdm=True,
):
    """Get genes for given descriptions using asyncio.

    Allows experiments with role prompt, different models,
    use of confidence and model reasoning. Used also for ensembling.

    Args:
       aclient: async OpenAI client
       descr: list of pathway/process descriptions
       model: OpenAI model string
       prompt_type: "basic", standard prompt, "reason",
                     add reasoning per gene, "conf" add confidence
                     per gene
       use_sysmsg: apply system message to model input
       n_retry: number of times to retry
       use_tqdm: true/false show tqdm progress bar
    Returns:
      list of a list of dicts with
      genes, unique genes in gene set
      parsed_genes, all genes parsed, including any duplicates
      reason, if requested reason gene was included in gene set, one for each
              gene in parse_genes, otherwise an empty string
      in_toks, input token count, out_toks, output count
      ntries, number of tries to get a gene set
    """
    prompt_file = "genes_concise.txt"

    # If requested, use prompts that require reasoning or confidence.
    if prompt_type == "reason":
        prompt_file = "genes_concise_reason.txt"
    if prompt_type == "conf":
        prompt_file = "genes_concise_conf.txt"

    with resources.open_text("llm2geneset.prompts", prompt_file) as file:
        prompt = file.read()

    encoding = tiktoken.encoding_for_model(model)

    sys_msg = "You are an expert in cellular and molecular biology."

    prompts = [prompt.format(descr=d) for d in descr]

    # TODO: Make this rate limit a parameter.
    rate_limiter = StrictLimiter(0.95 * 10000.0 / 60.0)

    async def complete(p):
        await rate_limiter.wait()
        in_toks = 0
        out_toks = 0
        for attempt in range(n_retry):
            # Count input tokens.
            in_toks += len(encoding.encode(p))
            # Prepend sys message if requested.
            messages = [{"role": "user", "content": p}]
            if use_sysmsg:
                messages = [{"role": "system", "content": sys_msg}] + messages
                in_toks += len(encoding.encode(sys_msg))
            # LLM
            r = await aclient.chat.completions.create(model=model, messages=messages)
            resp = r.choices[0].message.content
            # Count output tokens.
            out_toks += len(encoding.encode(resp))
            # Extract gene names.
            genes = []
            try:
                last_code = extract_last_code_block(resp)
                json_parsed = json_repair.loads(last_code)
                # Address issue where sometimes other types are parsed out.
                json_parsed = [g for g in json_parsed if isinstance(g["gene"], str)]
                # Parse out gene, reason, and confidence.
                genes = [g["gene"] for g in json_parsed]
                reason = ["" for g in json_parsed]
                if prompt_type == "reason":
                    reason = [g["reason"] for g in json_parsed]
                conf = ["" for g in json_parsed]
                if prompt_type == "conf":
                    conf = [g["confidence"] for g in json_parsed]
                return {
                    "parsed_genes": genes,
                    "reason": reason,
                    "conf": conf,
                    "in_toks": in_toks,
                    "out_toks": out_toks,
                    "ntries": attempt + 1,
                }
            except Exception as e:
                if attempt == n_retry - 1:
                    raise RuntimeError("Retries exceeded.") from e

    # Run completions asynchronously. Display progress bar if requested.
    if use_tqdm:
        res = await tqdm.asyncio.tqdm.gather(*(complete(p) for p in prompts))
    else:
        res = await asyncio.gather(*(complete(p) for p in prompts))
    return res


async def get_genes_context(
    aclient, context: List[str], descr: List[str], model="gpt-4o", n_retry=3
):
    """Get genes using given context.

    Generates a gene set given some context text.

     Args:
       aclient: async OpenAI client
       descr: list of natural language descriptions of a gene set
       context: list of textual context, a string
       model: OpenAI model name
       n_retry: number of retries per

    """
    prompt_file = "genes_concise_context.txt"

    with resources.open_text("llm2geneset.prompts", prompt_file) as file:
        prompt = file.read()

    prompt_file_orig = "genes_concise.txt"
    with resources.open_text("llm2geneset.prompts", prompt_file_orig) as file:
        prompt_orig = file.read()

    encoding = tiktoken.encoding_for_model(model)

    # If no context, use original prompt.
    prompts = []
    for c, d in zip(context, descr):
        if len(c) == 0:
            prompts.append(prompt_orig.format(descr=d))
        else:
            prompts.append(prompt.format(context=c, descr=d))

    # TODO: Make this rate limit a parameter.
    rate_limiter = StrictLimiter(0.95 * 10000.0 / 60.0)

    async def complete(p):
        await rate_limiter.wait()
        in_toks = 0
        out_toks = 0
        for attempt in range(n_retry):
            # Count input tokens.
            in_toks += len(encoding.encode(p))
            # Prepend sys message if requested.
            messages = [{"role": "user", "content": p}]
            # LLM
            r = await aclient.chat.completions.create(model=model, messages=messages)
            resp = r.choices[0].message.content
            # Count output tokens.
            out_toks += len(encoding.encode(resp))
            # Extract gene names.
            genes = []
            try:
                last_code = extract_last_code_block(resp)
                json_parsed = json_repair.loads(last_code)
                # Address issue where sometimes other types are parsed out.
                json_parsed = [g for g in json_parsed if isinstance(g["gene"], str)]
                # Parse out gene, reason, and confidence.
                genes = [g["gene"] for g in json_parsed]
                reason = ["" for g in json_parsed]
                conf = ["" for g in json_parsed]
                return {
                    "parsed_genes": genes,
                    "reason": reason,
                    "conf": conf,
                    "in_toks": in_toks,
                    "out_toks": out_toks,
                    "ntries": attempt + 1,
                }
            except Exception as e:
                print("retrying")
                print(e)
                print(p)
                print(resp)
                if attempt == n_retry - 1:
                    raise RuntimeError("Retries exceeded.") from e

    # Run completions asynchronously.
    res = await tqdm.asyncio.tqdm.gather(*(complete(p) for p in prompts))
    return res


def filter_items_by_threshold(list_of_lists, threshold):
    """Find repeated items in a list of lists.

    Args:
       list_of_lists: list of lists
       threshold: integer threshold for element
    Returns:
       list of unique elements that occur in threshold
       number of lists.
    """
    item_count = {}

    # Count the number of lists each item appears in
    for sublist in list_of_lists:
        unique_items = set(sublist)
        for item in unique_items:
            if item in item_count:
                item_count[item] += 1
            else:
                item_count[item] = 1

    # Filter items based on the threshold
    res_filtered = [item for item, count in item_count.items() if count >= threshold]
    result = sorted(res_filtered)

    return result


def ensemble_genes(descr, gen_genes, thresh):
    """Ensemble gene sets.

    Uses multiple generations of get_genes() to create gene sets.

    Args:
       descr: list of gene set descriptions
       gen_genes: output of a list of dicts from get_genes
       thresh: integer for how many generations a gene needs to
               appear
    Returns:
       Returns a dict with common genes in
       "genes", tokens are summed along with number of tries
       needed to generate the gene set.
    """
    ensembl_genes = []
    for idx in range(len(descr)):
        gene_lists = []
        in_toks = 0
        out_toks = 0
        ntries = 0
        for e in range(len(gen_genes)):
            gene_lists.append(gen_genes[e][idx]["parsed_genes"])
            in_toks += gen_genes[e][idx]["in_toks"]
            out_toks += gen_genes[e][idx]["out_toks"]
            ntries += gen_genes[e][idx]["ntries"]

        thresh_genes = filter_items_by_threshold(gene_lists, thresh)
        blank_list = ["" for g in thresh_genes]
        x = {
            "genes": thresh_genes,
            "parsed_genes": thresh_genes,
            "reason": blank_list,
            "conf": blank_list,
            "in_toks": in_toks,
            "out_toks": out_toks,
            "ntries": ntries,
        }
        ensembl_genes.append(x)

    return ensembl_genes


def sel_conf(descr, gen_genes, conf_vals):
    """
    Select genes based on given confidence values.

    Args:
       descr: list of gene set descriptions
       gen_genes: output of a list of dicts from get_genes
       thresh: integer for how many generations a gene needs to
               appear
    Returns:
       Returns a dict with common genes in
       "genes", tokens are summed along with number of tries
       needed to generate the gene set.
    """
    conf_vals = set(conf_vals)
    conf_genes = []
    for idx in range(len(descr)):
        genes = gen_genes[idx]["parsed_genes"]
        conf = gen_genes[idx]["conf"]
        reason = gen_genes[idx]["reason"]

        genes_sel = []
        conf_sel = []
        reason_sel = []
        for g in range(len(genes)):
            if conf[g] in conf_vals:
                genes_sel.append(genes[g])
                conf_sel.append(conf[g])
                reason_sel.append(reason[g])

        x = {
            "parsed_genes": genes_sel,
            "reason": reason_sel,
            "conf": conf_sel,
            "in_toks": gen_genes[idx]["in_toks"],
            "out_toks": gen_genes[idx]["out_toks"],
            "ntries": gen_genes[idx]["ntries"],
        }
        conf_genes.append(x)

    return conf_genes


async def gsai(aclient, protein_lists: List[List[str]], model="gpt-4o", n_retry=3):
    """Run GSAI from Ideker Lab.

    Uses the prompt from: https://idekerlab.ucsd.edu/gsai/ to summarize genes and
    uncover their function, also provide confidence.

    Args:
       aclient: asynchronous OpenAI client
       protein_lists: list of a list of genes, gene sets to
                 assign function
       model: OpenAI model string
       n_retry: number of retries to get valid parsed output
    """
    prompt_file = "gsai_prompt.txt"
    with resources.open_text("llm2geneset.prompts", prompt_file) as file:
        prompt = file.read()

    prompts = [prompt.format(proteins=", ".join(p)) for p in protein_lists]
    rate_limiter = StrictLimiter(0.95 * 10000.0 / 60.0)
    encoding = tiktoken.encoding_for_model(model)
    sys_msg = "You are an efficient and insightful assistant to a molecular biologist."

    # TODO: Make this rate limit a parameter.
    rate_limiter = StrictLimiter(0.95 * 10000.0 / 60.0)

    def parse_name(text):
        pattern = r"Name:\s*(.+?)\n"
        nmatch = re.search(pattern, text)
        if nmatch:
            return nmatch.group(1)
        else:
            return None

    def parse_conf(text):
        pattern = r"LLM self-assessed confidence:\s*([\d\.]+)"
        cmatch = re.search(pattern, text)
        if cmatch:
            return float(cmatch.group(1))
        else:
            return None

    def parse_list(input_text):
        list_items = re.findall(
            r"^\d+\.\s.*?(?=\n|\Z)", input_text, re.DOTALL | re.MULTILINE
        )
        text_after_list = [re.sub(r"^\d+\.\s", "", item) for item in list_items]
        text_after_list = [item.strip() for item in text_after_list]
        return list_items

    async def complete(p):
        await rate_limiter.wait()
        in_toks = 0
        out_toks = 0
        for attempt in range(n_retry):
            # Count input tokens.
            in_toks += len(encoding.encode(sys_msg))
            in_toks += len(encoding.encode(p))
            # Generate message.
            messages = [
                {"role": "system", "content": sys_msg},
                {"role": "user", "content": p},
            ]
            # LLM
            r = await aclient.chat.completions.create(model=model, messages=messages)
            resp = r.choices[0].message.content
            # Count output tokens.
            out_toks += len(encoding.encode(resp))
            try:
                name = parse_name(resp)
                conf = parse_conf(resp)
                annot = parse_list(resp)
                if name is None:
                    raise ValueError("name is none")
                if conf is None:
                    raise ValueError("conf is none")
                if annot is None:
                    raise ValueError("annot is none")
                return {
                    "name": name,
                    "conf": conf,
                    "annot": annot,
                    "in_toks": in_toks,
                    "out_toks": out_toks,
                    "ntries": attempt + 1,
                }
            except Exception as e:
                print("retrying")
                print(e)
                print(p)
                print(resp)
                if attempt == n_retry - 1:
                    raise RuntimeError("Retries exceeded.") from e

    # Run completions asynchronously.
    res = await tqdm.asyncio.tqdm.gather(*(complete(p) for p in prompts))
    return res


async def bp_from_genes(
    aclient, model, genes: List[str], n_pathways=5, context="", n_retry=3
):
    """Propose a list of biological processes from a set of genes.

    Args:
       aclient: asynchronous OpenAI client
       model: OpenAI model
       genes: list of genes to use to propose
       context: string describing context to consider when proposing pathways.
       n_pathways: number of pathways to propose
       n_retry: number of retries to get correctly parsing output
    Returns:
       List of processes and pathways proposed based on input
       genes.

    """
    # Generate message.
    prompt_file = "pathways_from_genes.txt"
    if context != "":
        prompt_file = "pathways_from_genes_context.txt"

    with resources.open_text("llm2geneset.prompts", prompt_file) as file:
        prompt = file.read()

    # Create the prompts by formatting the template
    if context == "":
        p = prompt.format(n_pathways=n_pathways, genes=",".join(genes))
    else:
        p = prompt.format(context=context, n_pathways=n_pathways, genes=",".join(genes))

    print(p)
    encoding = tiktoken.encoding_for_model(model)
    in_toks = 0
    out_toks = 0
    for attempt in range(n_retry):
        in_toks += len(encoding.encode(p))
        messages = [{"role": "user", "content": p}]
        r = await aclient.chat.completions.create(model=model, messages=messages)
        resp = r.choices[0].message.content
        out_toks += len(encoding.encode(resp))
        try:
            last_code = extract_last_code_block(resp)
            json_parsed = json_repair.loads(last_code)  # Use json.loads directly
            json_parsed = [path for path in json_parsed if isinstance(path["p"], str)]
            pathways = [path["p"] for path in json_parsed]
            if len(pathways) == 0:
                raise ValueError("No pathways returned.")
            return {"pathways": pathways, "in_toks": in_toks, "out_toks": out_toks}
        except Exception as e:
            print("retrying")
            print(p)
            print(e)
            print(resp)
            if attempt == n_retry - 1:
                raise RuntimeError("Retries exceeded.") from e


async def gs_proposal(
    aclient,
    protein_lists: List[List[str]],
    model="gpt-4o",
    context="",
    n_background=19846,
    n_pathways=5,
    n_retry=1,
):
    """Proposal-based approach to map from genes to function.

    Args:
       aclient: asynchronous OpenAI client
       protein_lists: list of a list of genes, gene sets to
                      assign function
       model: OpenAI model string
       n_background: number of genes in background set
       n_pathways: number of pathways to propose given a gene list
       n_retry: number of retries to get valid parsed output
    Returns:
      A dict with tot_in_toks (input) and tot_out_toks (output)
      tokens used. A pandas data frame with the hypergeometric
      overrepresentation results for each proposed gene set.
    """
    rate_limiter = StrictLimiter(0.95 * 10000.0 / 60.0)

    async def gse(genes):
        await rate_limiter.wait()
        # 1. Examine genes and propose possible pathways and processes.
        bio_process = await bp_from_genes(aclient, model, genes, n_pathways, context)

        # 2. Generate these gene sets without input genes as context.
        proposed = await get_genes(
            aclient, bio_process["pathways"], model=model, use_tqdm=False
        )

        # 3. Get over-representation p-values.
        tot_in_toks = bio_process["in_toks"]
        tot_out_toks = bio_process["out_toks"]
        output = []
        for idx in range(len(bio_process["pathways"])):
            llm_genes = list(set(proposed[idx]["parsed_genes"]))
            # Use hypergeometric to compute p-value.
            intersection = set(llm_genes).intersection(set(genes))
            p_val = hypergeom.sf(
                len(intersection) - 1,
                n_background - len(genes),
                len(genes),
                len(llm_genes),
            )
            tot_in_toks += proposed[idx]["in_toks"]
            tot_out_toks += proposed[idx]["out_toks"]

            generatio = float(len(intersection)) / len(set(genes))
            bgratio = float(len(set(llm_genes))) / n_background

            richFactor = None
            foldEnrich = None
            if len(llm_genes) > 0:
                richFactor = float(len(intersection)) / len(set(llm_genes))
                foldEnrich = generatio / bgratio

            output.append(
                {
                    "set_descr": bio_process["pathways"][idx],
                    "generatio": generatio,
                    "bgratio": bgratio,
                    "richFactor": richFactor,
                    "foldEnrich": foldEnrich,
                    "p_val": p_val,
                    "intersection": ",".join(list(intersection)),
                    "set_genes": ",".join(llm_genes),
                    "ngenes": len(set(genes)),
                    "nset": len(llm_genes),
                    "ninter": len(intersection),
                    "in_toks": proposed[idx]["in_toks"],
                    "out_toks": proposed[idx]["out_toks"],
                }
            )
        # Generate output, adjust p-values.
        df = pd.DataFrame(output)
        df.sort_values("p_val", inplace=True)
        _, p_adj, _, _ = multipletests(df["p_val"], method="fdr_bh")
        df["p_adj"] = p_adj
        loc = df.columns.get_loc("p_val") + 1
        new_columns = df.columns.tolist()
        new_columns.insert(loc, new_columns.pop(new_columns.index("p_adj")))
        df = df[new_columns]
        return {
            "tot_in_toks": tot_in_toks,
            "tot_out_toks": tot_out_toks,
            "ora_results": df,
        }

    res = await tqdm.asyncio.tqdm.gather(*(gse(p) for p in protein_lists))
    return res


def simple_ora(genes: List[str], set_descr, gene_sets, n_background=19846):
    """
    Run simple overrepresentation analysis on a set of genes.

    Args:
       genes: genes on which to perform overreprsentation analysis
       set_descr: list of gene set descriptions
       gene_sets: list of a list of genes
       n_background: size of background gene set
    """
    output = []
    for idx in range(len(set_descr)):
        set_genes = list(set(gene_sets[idx]))
        # Use hypergeometric to compute p-value.
        intersection = set(set_genes).intersection(set(genes))
        p_val = hypergeom.sf(
            len(intersection) - 1,
            n_background - len(genes),
            len(genes),
            len(set_genes),
        )
        generatio = float(len(intersection)) / len(set(genes))
        bgratio = float(len(set(set_genes))) / n_background

        richFactor = None
        foldEnrich = None
        if len(set_genes) > 0:
            richFactor = float(len(intersection)) / len(set(set_genes))
            foldEnrich = generatio / bgratio

        output.append(
            {
                "set_descr": set_descr[idx],
                "generatio": generatio,
                "bgratio": bgratio,
                "richFactor": richFactor,
                "foldEnrich": foldEnrich,
                "p_val": p_val,
                "intersection": ",".join(list(intersection)),
                "set_genes": ",".join(set_genes),
                "ngenes": len(set(genes)),
                "nset": len(set_genes),
                "ninter": len(intersection),
            }
        )

    # Generate output, adjust p-values.
    df = pd.DataFrame(output)
    df.sort_values("p_val", inplace=True)
    _, p_adj, _, _ = multipletests(df["p_val"], method="fdr_bh")
    df["p_adj"] = p_adj
    loc = df.columns.get_loc("p_val") + 1
    new_columns = df.columns.tolist()
    new_columns.insert(loc, new_columns.pop(new_columns.index("p_adj")))
    df = df[new_columns]
    return df
