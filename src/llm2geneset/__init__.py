import sys

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

from .llm2geneset import (
    ensemble_genes,
    get_embeddings,
    get_genes,
    get_genes_context,
    gs_proposal,
    gsai,
    read_gmt,
    sel_conf,
)

__all__ = [
    "read_gmt",
    "get_genes",
    "ensemble_genes",
    "sel_conf",
    "get_genes_context",
    "get_embeddings",
    "gsai",
    "gs_proposal",
]


from .eutils import efetch_pubmed_async, efetch_pubmed_sync, esearch_async, esearch_sync

__all__ += [
    "esearch_async",
    "esearch_sync",
    "efetch_pubmed_async",
    "efetch_pubmed_sync",
]
