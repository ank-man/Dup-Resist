"""Dup-Resist: utilities for identifying gene conversion among duplicated genes."""

from .data_models import GenePair, CollinearityBlock, GeneConversionEvent
from .mcscanx_parser import (
    parse_gene_types,
    parse_collinearity,
    parse_kaks_table,
)
from .gene_conversion import GeneConversionDetector
from .singletonfx import (
    SingletonFXPipeline,
    build_singletonfx_parser,
    run_singletonfx,
)

__all__ = [
    "GenePair",
    "CollinearityBlock",
    "GeneConversionEvent",
    "GeneConversionDetector",
    "parse_gene_types",
    "parse_collinearity",
    "parse_kaks_table",
    "SingletonFXPipeline",
    "build_singletonfx_parser",
    "run_singletonfx",
]
