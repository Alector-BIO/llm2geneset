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

This project enables generation of gene sets using natural language descriptions
of biological pathways and processes.

## Usage

Please read
[OpenAI's best practicies for API key safety](https://help.openai.com/en/articles/5112595-best-practices-for-api-key-safety)
and make sure you have your OpenAI API key setup as
environment variables e.g.:

```bash
export OPENAI_API_ORG="org-XXXX"
export OPENAI_API_KEY="XXXXX"
```

It is also recommended you obtain an [NCBI API Key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
and add it to your environment.

```bash
export NCBI_API_KEY="x"
```

```python
import llm2geneset

# New interface is TBD.
```

## Installation using micromamba

Create an environment using micromamba and work on this
package in editable mode. Here's how to install micromamba:

```bash
cd ~
# Linux x86_64
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
# M1 Mac Silicon
curl -Ls https://micro.mamba.pm/api/micromamba/osx-arm64/latest | tar -xvj bin/micromamba
```

In your .bashrc add `bin` to your path:

```bash
export PATH="$HOME/bin:${PATH}"
```

Then run the following to setup a location for all of the micromamba environments.

```bash
micromamba shell init -s bash -p ~/micromamba
```

**NOTE: The rust compiler is only for GSEApy which is in the environment. We
may or may not keep this dependency.**

Install the main environment. Note that you need `rustc` the
rust compiler installed. You can install `rustc` by running
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
Then, add `. "$HOME/.cargo/env"` to your .bashrc.

```bash
micromamba env create -f llm2geneset.yml
micromamba activate llm2geneset
```

Install the package in editable mode.

```bash
cd llm2geneset
pip install -e .
```


<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
