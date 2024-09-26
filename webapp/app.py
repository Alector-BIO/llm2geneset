"""The llm2geneset streamlit app."""

import asyncio
import os

import nest_asyncio
import openai

# import pandas as pd
import streamlit as st

import llm2geneset

nest_asyncio.apply()


async def async_function(aclient, genes, model, context):
    res = await llm2geneset.gs_proposal(
        aclient, genes, model=model, context=context, n_pathways=100
    )
    return res


with st.sidebar:
    default_key = os.getenv("OPENAI_API_KEY")
    if default_key:
        openai_api_key = st.text_input("OpenAI API Key", value=default_key)
    else:
        openai_api_key = st.text_input("OpenAI API Key", value=default_key)
    "[Get an OpenAI API key](https://platform.openai.com)"
    "[View source](https://github.com/Alector-Biotech/llm2geneset)"

st.title("llm2geneset")
st.caption("Demonstration app for llm2geneset.")

model_input = st.text_area("OpenAI Model", value="gpt-4o-2024-05-13")

# Text area to input experimental context
context_input = st.text_area(
    "Experimental Context", placeholder="Enter experimental context here..."
)

# Text area to input genes
genes_input = st.text_area("Gene List", placeholder="Enter genes here...")

# Button to trigger the action
if st.button("Go"):
    # Split the input by commas or newlines
    gene_list = genes_input.replace(",", "\n").splitlines()

    # Process the genes
    genes = [gene.strip() for gene in gene_list if gene.strip()]
    model = model_input.strip()
    context = context_input.strip()

    # Display the processed genes
    st.write("Context:")
    st.write(context)
    st.write("Genes:")
    st.write(genes)

    with st.spinner("Running llm2geneset please wait..."):
        # Use asyncio to run the asynchronous function
        aclient = openai.AsyncClient()
        loop = asyncio.get_event_loop()
        afun = async_function(aclient, genes, model, context)
        res = loop.run_until_complete(afun)
        in_toks = res["tot_in_toks"]
        out_toks = res["tot_out_toks"]
        st.write(f"input tokens {in_toks}, output tokens {out_toks}")
        st.dataframe(res["ora_results"])
    st.success("Done!")
