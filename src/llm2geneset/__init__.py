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
    bp_from_genes,
    ensemble_genes,
    get_embeddings,
    get_genes,
    get_genes_bench,
    gs_proposal,
    gs_proposal_bench,
    gsai_bench,
    read_gmt,
    sel_conf,
    simple_ora,
)

__all__ = [
    "bp_from_genes",
    "read_gmt",
    "get_embeddings",
    "get_genes_bench",
    "ensemble_genes",
    "get_genes",
    "sel_conf",
    "gsai_bench",
    "gs_proposal_bench",
    "gs_proposal",
    "simple_ora",
]
