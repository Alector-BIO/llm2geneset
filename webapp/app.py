"""The llm2geneset streamlit app."""

import asyncio
import os

import nest_asyncio
import openai

# import pandas as pd
import streamlit as st

import llm2geneset

nest_asyncio.apply()


async def async_function(
    aclient, genes, model, context, n_pathways, n_background, seed
):
    res = await llm2geneset.gs_proposal(
        aclient,
        genes,
        model=model,
        context=context,
        n_pathways=n_pathways,
        n_background=n_background,
        #bgd_genes = bgd_genes,
        seed=seed,
    )
    return res


with st.sidebar:
    default_key = os.getenv("OPENAI_API_KEY")
    if default_key:
        openai_api_key = st.text_input("OpenAI API Key", value=default_key, type="password")
    else:
        openai_api_key = st.text_input("OpenAI API Key", value="", type="password")
    "[Get an OpenAI API key](https://platform.openai.com)"
    "[View source](https://github.com/Alector-Biotech/llm2geneset)"

st.title("llm2geneset")
st.caption("Demonstration app for llm2geneset.")

model_input = st.text_area("OpenAI Model", value="gpt-4o-2024-05-13")

genes_input = st.text_area("Gene List", placeholder="Enter genes here...")

context_input = st.text_area(
    "Experimental Context", placeholder="Enter experimental context here..."
)

num_gene_sets = st.number_input("# of gene sets", value=100, min_value=0)

n_background = st.number_input("# of background genes ", value=19846, min_value=0)

bgd_genes = None

seed = st.number_input("seed", value=3272995, min_value=0)

# Button to trigger the action
if st.button("Go"):
    # Split the input by commas or newlines
    gene_list = genes_input.replace(",", "\n").splitlines()

    # Process the genes
    genes = [gene.strip() for gene in gene_list if gene.strip()]
    model = model_input.strip()
    context = context_input.strip()

    # Display the processed genes
    if context != "":
        st.write("Context:")
        st.write(context)

    st.write("Genes:")
    st.write(genes)

    with st.spinner("Running llm2geneset please wait..."):
        # Use asyncio to run the asynchronous function
        aclient = openai.AsyncClient(api_key=openai_api_key)
        loop = asyncio.get_event_loop()
        afun = async_function(
            aclient, genes, model, context, num_gene_sets, n_background, seed
        )
        res = loop.run_until_complete(afun)
        in_toks = res["tot_in_toks"]
        out_toks = res["tot_out_toks"]
        st.write(f"input tokens {in_toks}, output tokens {out_toks}")
        st.dataframe(res["ora_results"])
    st.success("Done!")
