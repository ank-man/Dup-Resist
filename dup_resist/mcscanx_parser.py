"""Parsers for MCScanX output files."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .data_models import CollinearityBlock, GenePair


def parse_gene_types(path: Path) -> Dict[str, str]:
    """Parse the ``*.gene_type`` file produced by ``duplicate_gene_classifier``.

    Parameters
    ----------
    path:
        Path to the ``*.gene_type`` file.

    Returns
    -------
    Dict[str, str]
        Mapping from gene identifier to duplicate class (e.g. ``WGD``).
    """

    gene_types: Dict[str, str] = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            gene_types[parts[0]] = parts[1]
    return gene_types


def parse_collinearity(path: Path) -> List[CollinearityBlock]:
    """Parse an ``*.collinearity`` file produced by MCScanX."""

    blocks: List[CollinearityBlock] = []
    current_block: Optional[CollinearityBlock] = None
    pair_index = 0

    with Path(path).open() as handle:
        for raw in handle:
            line = raw.rstrip()
            if not line:
                continue
            if line.startswith("#"):  # comment or alignment header
                if line.startswith("## Alignment"):
                    # Example: "## Alignment 1: score=134.0 e_value=1e-30 N=14" etc
                    parts = line.split()
                    block_id = parts[2].rstrip(":")
                    score = _parse_key_value(line, "score")
                    e_value = _parse_key_value(line, "e_value")
                    chrom_info = _extract_chromosomes(line)
                    current_block = CollinearityBlock(
                        block_id=block_id,
                        score=score,
                        e_value=e_value,
                        chromosome1=chrom_info[0],
                        chromosome2=chrom_info[1],
                    )
                    blocks.append(current_block)
                    pair_index = 0
                continue

            if current_block is None:
                continue

            pair_index += 1
            parts = line.split()
            if len(parts) < 5:
                # Usually contains match score, orientation, gene ids, etc.
                continue

            try:
                score = float(parts[0])
            except ValueError:
                score = None
            orientation = parts[1]
            anchor1 = parts[2]
            anchor2 = parts[3]
            try:
                e_value = float(parts[4])
            except ValueError:
                e_value = None

            pair = GenePair(
                anchor1=anchor1,
                anchor2=anchor2,
                block_id=current_block.block_id,
                pair_index=pair_index,
                score=score,
                e_value=e_value,
                orientation=orientation,
            )
            current_block.add_pair(pair)

    return blocks


def _extract_chromosomes(line: str) -> Tuple[Optional[str], Optional[str]]:
    chrom1: Optional[str] = None
    chrom2: Optional[str] = None
    if "between" in line:
        # Example: "between Chr1 and Chr3"
        fragments = line.split("between", 1)[1].strip()
        parts = fragments.split()
        if len(parts) >= 3:
            chrom1 = parts[0]
            chrom2 = parts[2]
    return chrom1, chrom2


def _parse_key_value(line: str, key: str) -> Optional[float]:
    token = f"{key}="
    if token not in line:
        return None
    value = line.split(token, 1)[1].split()[0]
    try:
        return float(value)
    except ValueError:
        return None


def parse_kaks_table(path: Path) -> Dict[Tuple[str, str], Tuple[Optional[float], Optional[float]]]:
    """Parse the ``*.kaks.tsv`` table emitted by ``add_ka_and_ks_to_collinearity.pl``."""

    kaks: Dict[Tuple[str, str], Tuple[Optional[float], Optional[float]]] = {}
    with Path(path).open() as handle:
        header: Optional[List[str]] = None
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if header is None:
                header = parts
                continue

            row = dict(zip(header, parts))
            gene1 = row.get("gene1") or row.get("anchor1") or parts[0]
            gene2 = row.get("gene2") or row.get("anchor2") or parts[1]
            ka = _safe_float(row.get("Ka") or row.get("ka"))
            ks = _safe_float(row.get("Ks") or row.get("ks"))

            key = (gene1, gene2)
            kaks[key] = (ka, ks)
            # Store both orientations for convenience
            kaks[(gene2, gene1)] = (ka, ks)
    return kaks


def _safe_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None
