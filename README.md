<!-- These are examples of badges you might want to add to your README:
     please update the URLs accordingly

[![Built Status](https://api.cirrus-ci.com/github/<USER>/llm2geneset.svg?branch=main)](https://cirrus-ci.com/github/<USER>/llm2geneset)
[![ReadTheDocs](https://readthedocs.org/projects/llm2geneset/badge/?version=latest)](https://llm2geneset.readthedocs.io/en/stable/)
[![Coveralls](https://img.shields.io/coveralls/github/<USER>/llm2geneset/main.svg)](https://coveralls.io/r/<USER>/llm2geneset)
[![PyPI-Server](https://img.shields.io/pypi/v/llm2geneset.svg)](https://pypi.org/project/llm2geneset/)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/llm2geneset.svg)](https://anaconda.org/conda-forge/llm2geneset)
[![Monthly Downloads](https://pepy.tech/badge/llm2geneset/month)](https://pepy.tech/project/llm2geneset)
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter)](https://twitter.com/llm2geneset)
-->

[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)

# llm2geneset

> llm2geneset

This project combines LLMs, gene set generation, and overrepresentation analysis
to power analyis of RNA-seq, scRNA-seq, and proteomics data sets.

llm2geneset is similar to popular tools such as [enrichr](https://maayanlab.cloud/Enrichr/),
[cluststerProfiler](https://guangchuangyu.github.io/software/clusterProfiler/), or
[DAVID](https://davidbioinformatics.nih.gov/), but is uses LLMs to propose gene set descriptions
and gene sets themselves. The generated gene sets can also be used in tools such as
[fGSEA](https://bioconductor.org/packages/release/bioc/html/fgsea.html) and
[GSEApy](https://gseapy.readthedocs.io/en/latest/introduction.html).

If you have an OpenAI API key, you can try out the web application
at [https://llm2geneset.streamlit.app](https://llm2geneset.streamlit.app).

## Link to the Paper
[Enhancing Gene Set Overrepresentation Analysis with Large Language Models](https://doi.org/10.1093/bioadv/vbaf054)

Jiqing Zhu, Rebecca Y. Wang, Xiaoting Wang, Ricardo Azevedo, Alex Moreno, Julia A. Kuhn,  Zia Khan


## OpenAI API Key Setup

Please read
[OpenAI's best practicies for API key safety](https://help.openai.com/en/articles/5112595-best-practices-for-api-key-safety)
and make sure you have your OpenAI API key setup as
environment variables e.g.:

```bash
export OPENAI_API_ORG="org-XXXX"
export OPENAI_API_KEY="XXXXX"
```

## Installation using pixi

An environment to run `llm2geneset` can be configured using 
[pixi](https://prefix.dev)

Once pixi is installed, run `pixi shell` in the llm2geneset directory 

```bash
cd llm2geneset
pixi shell
```

## Installing as a pypi package

You can also install the package
using `pip` as it is available on pypi [https://pypi.org/project/llm2geneset](https://pypi.org/project/llm2geneset)

```bash
pip install llm2geneset
```

## Usage

You can use the package in a script as follows after 
running `pixi shell`.

```python
import openai
import llm2geneset
import asyncio

async def main():
    aclient = openai.AsyncClient()
    genes = await llm2geneset.get_genes(aclient, "Antigen presentation")
    print(','.join(genes['parsed_genes']))
    res = await llm2geneset.gs_proposal(
        aclient, genes, n_pathways=5,
        n_background=19846)
    res = await llm2geneset.gs_proposal(aclient, genes['parsed_genes'])
    print(res['ora_results'])

if __name__ == "__main__":
    asyncio.run(main())

# Output:
# HLA-A,HLA-B,HLA-C,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-DRB4,...
# set_descr  generatio   bgratio  richFactor  foldEnrich
# 1  Antigen processing and presentation via MHC cl...   0.500000  0.001209    0.625000  413.458333
# ...
```


See also [notebooks/simple_example.ipynb](notebooks/simple_example.ipynb) for an example
of how to use it within a jupyter notebook.

## Webapp Interface

Streamlit is included in the llm2geneset environment. You
can run the webapp interface as follows.

```bash
pixi shell
streamlit run webapp/app.py
```

## pypi deployment

Bump version in [setup.cfg](setup.cfg) and run the following according to the instructions
for [pyscaffold-markdown](https://github.com/pyscaffold/pyscaffoldext-markdown)
```bash
pixi shell
tox -e docs  # to build documentation
tox -e build -- --wheel # to build the package distribution w/o source
tox -e publish  # to test project uploads correctly in test.pypi.org
tox -e publish -- --repository pypi  # release package to PyPI
```



<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
