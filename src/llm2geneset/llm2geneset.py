"""llm2geneset: using LLMs to generate gene sets."""

import re
from importlib import resources

import json_repair
import tiktoken
import tqdm.asyncio
from asynciolimiter import Limiter
from sklearn.metrics.pairwise import cosine_similarity


def read_gmt(gmt_file: str):
    """
    Load a GMT file.

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


def get_embeddings(client, text_list, model="text-embedding-3-large"):
    """Get embeddings using OpenAI API, processing in batches of 2048.

    Args:
        text_list: lists of texts to embed
        model: embedding model
    Returns:
        List of embeddings.
    """

    def chunks(lst, n):
        """Yield successive n-sized chunks from lst (sliding window)."""
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    # Lowercase all text entries
    text_list = [text.lower() for text in text_list]

    # Process in batches of 2048
    all_embeddings = []
    for batch in chunks(text_list, 2048):  # for each sliding window
        # we get the response, which includes the embedding vector
        response = client.embeddings.create(model=model, input=batch).data
        # Extract embeddings for the current batch (sliding window) and append to
        # the all_embeddings list
        batch_embeddings = [resp.embedding for resp in response]
        # we append the embeddings for this sliding window to all embeddings
        all_embeddings.extend(batch_embeddings)
    # for the text_list, convert sliding windows of text into sliding
    # windows of embeddings
    return all_embeddings


def get_csim(
    client,
    lib,
    sample1,
    sample2,
    spre="biological function of ",
    gpre="genes involved in ",
):
    """Get cosine similarity between gene set description and sample (w/ prefix)."""
    if not spre.endswith(" "):
        spre += " "
    if not gpre.endswith(" "):
        gpre += " "

    gene_sets = list(lib.keys())
    gene_sets_descriptive = [gpre + i for i in gene_sets]

    if sample2 is not None:
        merged = gene_sets_descriptive + [spre + sample1] + [spre + sample2]
        emb = get_embeddings(client, merged)
        # compute similarity between geneset list, NOTE
        csim1 = cosine_similarity(emb[: len(gene_sets)], [emb[len(gene_sets)]])
        csim2 = cosine_similarity(emb[: len(gene_sets)], [emb[len(gene_sets) + 1]])

        csim1 = dict(zip(gene_sets, csim1.squeeze()))
        csim2 = dict(zip(gene_sets, csim2.squeeze()))

        return (csim1, csim2)
    else:
        merged = gene_sets_descriptive + [sample1]
        emb = get_embeddings(client, merged)
        csim1 = cosine_similarity(emb[: len(gene_sets)], [emb[len(gene_sets)]])
        csim1 = dict(zip(gene_sets, csim1.squeeze()))
        return (csim1, None)


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
    aclient, descr, model="gpt-4o", prompt_type="basic", use_sysmsg=False, n_retry=3
):
    """Get genes for given descriptions using asyncio.

    Args:
       aclient: async OpenAI client
       descr: list of pathway/process descriptions
       model: OpenAI model string
       prompt_type: "basic", standard prompt, "reason",
                     add reasoning per gene, "conf" add confidence
                     per gene
       use_sysmsg: apply system message to model input
       n_retry: number of times to retry
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
    rate_limiter = Limiter(0.95 * 10000.0 / 60.0)

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
    result = [item for item, count in item_count.items() if count >= threshold]

    return result


def ensemble_genes(descr, gen_genes, thresh):
    """Ensemble gene sets.

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
        genes = gen_genes[idx]["parse_genes"]
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
