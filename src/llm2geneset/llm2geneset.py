"""llm2geneset: using LLMs to generate gene sets."""

import re
from importlib import resources

import json_repair
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
    """Extract last code block in output."""
    # Regular expression to find code blocks enclosed in triple backticks
    pattern = r"```(?:[\w]*\n)?([\s\S]*?)```"
    code_blocks = re.findall(pattern, markdown_text)
    # Return the last code block, if any
    return code_blocks[-1] if code_blocks else None


async def get_genes(aclient, descr, modelg="gpt-4-turbo", species="human", n_retry=3):
    """Get genes for given descriptions using asyncio."""
    # TODO: Need to add the ability to count tokens.

    with resources.open_text("llm2geneset.prompts", "genes_concise.txt") as file:
        prompt = file.read()

    sys_msg = "You are a skilled assistant to a molecular biologist."

    prompts = [prompt.format(descr=d) for d in descr]

    # TODO: Make this rate limit a parameter.
    rate_limiter = Limiter(0.9 * 10000.0 / 60.0)

    async def complete(p):
        await rate_limiter.wait()
        for attempt in range(n_retry):
            r = await aclient.chat.completions.create(
                model=modelg,
                messages=[
                    {"role": "system", "content": sys_msg},
                    {"role": "user", "content": p},
                ],
            )
            resp = r.choices[0].message.content
            genes = []
            try:
                genes = json_repair.loads(extract_last_code_block(resp))
                genes = [g["gene"] for g in genes]
                genes = list(set(genes))
                return genes
            except Exception:
                print("retrying")
                print(p)
                if attempt == n_retry - 1:
                    return None

    res = await tqdm.asyncio.tqdm.gather(*(complete(p) for p in prompts))
    return res
