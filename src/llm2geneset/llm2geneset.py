"""GSEAAI - Gene set enrichment analysis + LLMs."""

import itertools
import re
from importlib import resources
from typing import List

import json_repair
import tqdm.asyncio
from asynciolimiter import Limiter
from sklearn.metrics.pairwise import cosine_similarity


def extract_last_code_block(markdown_text):
    """Extract last code block in output."""
    # Regular expression to find code blocks enclosed in triple backticks
    pattern = r"```(?:[\w]*\n)?([\s\S]*?)```"
    code_blocks = re.findall(pattern, markdown_text)
    # Return the last code block, if any
    return code_blocks[-1] if code_blocks else None


def read_gmt(gmt_file: str, background_genes: List[str] = [], verbose=False):
    """
    Load a GMT file.

    Adapted from
    https://github.com/MaayanLab/blitzgsea/blob/87c9b62c1239df1a9005c8ee2d9d2c1eea9b55ca/blitzgsea/enrichr.py#L50

    Args:
       gmt_file:

    """
    with open(gmt_file, "r") as file:
        lines = file.readlines()
    library = {}
    background_set = {}
    if len(background_genes) > 1:
        background_genes = [x.upper() for x in background_genes]
        background_set = set(background_genes)
    for line in lines:
        sp = line.strip().split("\t")
        sp2 = [re.sub(",.*", "", value) for value in sp[2:]]
        sp2 = [x.upper() for x in sp2 if x]
        if len(background_genes) > 2:
            geneset = list(set(sp2).intersection(background_set))
            if len(geneset) > 0:
                library[sp[0]] = geneset
        else:
            if len(sp2) > 0:
                library[sp[0]] = sp2
    ugenes = list(set(list(itertools.chain.from_iterable(library.values()))))
    if verbose:
        print(
            "Library loaded. Library contains "
            + str(len(library))
            + " gene sets. "
            + str(len(ugenes))
            + " unique genes found."
        )
    return library


class GSEAAI:
    """GSEA AI."""

    sys_msg = "You are a skilled assistant to a molecular biologist."

    def __init__(self, client):
        """Initialize GSEAI.

        Args:
            client: OpenAI API client
        """
        self.client = client  # this is just the synchronous client

    def get_embeddings(self, text_list, model="text-embedding-3-large"):
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
            response = self.client.embeddings.create(model=model, input=batch).data
            # Extract embeddings for the current batch (sliding window) and append to
            # the all_embeddings list
            batch_embeddings = [resp.embedding for resp in response]
            # we append the embeddings for this sliding window to all embeddings
            all_embeddings.extend(batch_embeddings)
        # for the text_list, convert sliding windows of text into sliding
        # windows of embeddings
        return all_embeddings

    def get_csim(
        self,
        lib,
        sample1,
        sample2,
        spre="biological function of ",
        gpre="genes involved in ",
    ):
        """
        @param lib: gene set library, a dictionary, keys are names of pathways,
                    values are individual genes
        @param sample1: descriptions of tissue x, string
        @param sample2: descriptions of tissue y, string, may be None (if disease study)
        @param spre: string prefix, in order to get better vector representation for
                     sample 1/2, added this text, pinpointing that we want biological
                     focus
        @param gpre: geneset prefix, for each geneset in library, have
                     geneset names, similar to prefix to samples, added this text
        """
        """Get cosine similarity between gene set description and sample (w/ prefix)."""
        if not spre.endswith(" "):
            spre += " "
        if not gpre.endswith(" "):
            gpre += " "

        gene_sets = list(lib.keys())
        gene_sets_descriptive = [gpre + i for i in gene_sets]

        if sample2 is not None:
            merged = gene_sets_descriptive + [spre + sample1] + [spre + sample2]
            emb = self.get_embeddings(merged)
            # compute similarity between geneset list, NOTE
            csim1 = cosine_similarity(emb[: len(gene_sets)], [emb[len(gene_sets)]])
            csim2 = cosine_similarity(emb[: len(gene_sets)], [emb[len(gene_sets) + 1]])

            csim1 = dict(zip(gene_sets, csim1.squeeze()))
            csim2 = dict(zip(gene_sets, csim2.squeeze()))

            return (csim1, csim2)
        else:
            merged = gene_sets_descriptive + [sample1]
            emb = self.get_embeddings(merged)
            csim1 = cosine_similarity(emb[: len(gene_sets)], [emb[len(gene_sets)]])
            csim1 = dict(zip(gene_sets, csim1.squeeze()))
            return (csim1, None)

    async def get_genes(
        self, aclient, descr, modelg="gpt-4-turbo", species="human", n_retry=3
    ):
        """Get genes for given descriptions using asyncio."""
        if n_retry == 0:
            print("retry limit reached")
            return {}

        if len(descr) == 0:
            return {}

        with resources.open_text("gseaai.prompts", "genes_concise.txt") as file:
            prompt = file.read()

        if species == "yeast":
            # seems to work in principle (needs more testing)
            with resources.open_text(
                "gseaai.prompts", "genes_concise_yeast.txt"
            ) as file:
                prompt = file.read()

        prompts = [prompt.format(descr=d) for d in descr]

        rate_limiter = Limiter(0.9 * 10000.0 / 60.0)

        async def query(p):
            await rate_limiter.wait()
            return await aclient.chat.completions.create(
                model=modelg,
                messages=[
                    {"role": "system", "content": self.sys_msg},
                    {"role": "user", "content": p},
                ],
            )

        bres = await tqdm.asyncio.tqdm.gather(*(query(p) for p in prompts))

        # Extract pathway/bp and json responses.
        try:
            pathway_or_process_list = descr
            gene_response = [r.choices[0].message.content for r in bres]
        except Exception:
            pathway_or_process_list = []
            gene_response = []

        # Create library by converting gene list from json.
        lib = {}
        success_pathways = set()
        for p, genes_json in zip(pathway_or_process_list, gene_response):
            try:
                genes = json_repair.loads(extract_last_code_block(genes_json))
                genes = [g["gene"] for g in genes]
            except Exception as err:
                # Error in the json. Queue this up for retrying.
                print(err)
                print(p)
                print(genes_json)
            else:
                lib[p] = list(set(genes))
                if len(lib[p]) != len(genes):
                    print("warning: genes in " + p + " are not unique")
                success_pathways.add(p)

        # Retry based on failed pathways. This also handles the case where
        # an API request times out.
        failed_descr = {d for d in descr if d not in success_pathways}
        if len(failed_descr) > 0:
            print("failed pathways")
            print(failed_descr)

        failed = await self.get_genes(
            aclient, failed_descr, modelg, n_retry=n_retry - 1
        )

        return {**lib, **failed}
