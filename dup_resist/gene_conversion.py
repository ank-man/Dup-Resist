"""Algorithms for detecting gene conversion among duplicated genes."""

from __future__ import annotations

import statistics
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .data_models import CollinearityBlock, GeneConversionEvent, GenePair
from .mcscanx_parser import parse_collinearity, parse_gene_types, parse_kaks_table


class GeneConversionDetector:
    """Identify putative gene conversion events from MCScanX output."""

    def __init__(
        self,
        duplicate_class: str = "WGD",
        min_block_size: int = 4,
        z_threshold: float = 1.5,
        absolute_ks_threshold: Optional[float] = None,
        min_pairs_with_ks: int = 3,
    ) -> None:
        self.duplicate_class = duplicate_class
        self.min_block_size = min_block_size
        self.z_threshold = z_threshold
        self.absolute_ks_threshold = absolute_ks_threshold
        self.min_pairs_with_ks = min_pairs_with_ks

    def from_species_directory(
        self,
        species: str,
        species_dir: Path,
        output_dir: Path,
    ) -> List[GeneConversionEvent]:
        """Run detection for a single species directory."""

        gene_type_path = species_dir / f"{species}.gene_type"
        collinearity_path = species_dir / f"{species}.collinearity"
        kaks_path = species_dir / f"{species}.kaks.tsv"

        if not gene_type_path.exists():
            raise FileNotFoundError(gene_type_path)
        if not collinearity_path.exists():
            raise FileNotFoundError(collinearity_path)
        if not kaks_path.exists():
            raise FileNotFoundError(kaks_path)

        gene_types = parse_gene_types(gene_type_path)
        blocks = parse_collinearity(collinearity_path)
        kaks = parse_kaks_table(kaks_path)

        events = self.detect(species, blocks, gene_types, kaks)
        output_dir.mkdir(parents=True, exist_ok=True)
        self._write_events(output_dir / f"{species}_gene_conversion.tsv", events)
        return events

    def detect(
        self,
        species: str,
        blocks: Sequence[CollinearityBlock],
        gene_types: Dict[str, str],
        kaks: Dict[Tuple[str, str], Tuple[Optional[float], Optional[float]]],
    ) -> List[GeneConversionEvent]:
        """Detect gene conversion events across blocks."""

        events: List[GeneConversionEvent] = []
        duplicate_class = self.duplicate_class

        for block in blocks:
            wgd_pairs = [
                pair
                for pair in block
                if gene_types.get(pair.anchor1) == duplicate_class
                and gene_types.get(pair.anchor2) == duplicate_class
            ]
            if len(wgd_pairs) < self.min_block_size:
                continue

            decorated_pairs = self._attach_evolutionary_rates(wgd_pairs, kaks)
            pairs_with_ks = [p for p in decorated_pairs if p.ks is not None]
            if len(pairs_with_ks) < self.min_pairs_with_ks:
                continue

            ks_values = [p.ks for p in pairs_with_ks if p.ks is not None]
            median_ks = statistics.median(ks_values)
            stdev = statistics.pstdev(ks_values) if len(ks_values) > 1 else None

            for pair in pairs_with_ks:
                if pair.ks is None:
                    continue
                zscore = None
                if stdev and stdev > 0:
                    zscore = (pair.ks - median_ks) / stdev

                reason_parts: List[str] = []
                is_conversion = False

                if zscore is not None and zscore <= -self.z_threshold:
                    is_conversion = True
                    reason_parts.append(
                        f"Ks z-score {zscore:.2f} <= -{self.z_threshold:.2f}"
                    )

                if (
                    self.absolute_ks_threshold is not None
                    and pair.ks <= self.absolute_ks_threshold
                ):
                    is_conversion = True
                    reason_parts.append(
                        f"Ks {pair.ks:.3f} <= {self.absolute_ks_threshold:.3f}"
                    )

                if is_conversion:
                    events.append(
                        GeneConversionEvent(
                            species=species,
                            block_id=block.block_id,
                            pair_index=pair.pair_index,
                            gene1=pair.anchor1,
                            gene2=pair.anchor2,
                            ks=pair.ks,
                            ka=pair.ka,
                            block_median_ks=median_ks,
                            block_zscore=zscore,
                            detection_reason="; ".join(reason_parts) or "threshold",
                        )
                    )
        return events

    def _attach_evolutionary_rates(
        self,
        pairs: Iterable[GenePair],
        kaks: Dict[Tuple[str, str], Tuple[Optional[float], Optional[float]]],
    ) -> List[GenePair]:
        decorated: List[GenePair] = []
        for pair in pairs:
            ka, ks = kaks.get(pair.as_tuple(), (None, None))
            pair.ka = ka
            pair.ks = ks
            decorated.append(pair)
        return decorated

    def _write_events(self, path: Path, events: Sequence[GeneConversionEvent]) -> None:
        import csv

        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "species",
                    "block_id",
                    "pair_index",
                    "gene1",
                    "gene2",
                    "Ks",
                    "Ka",
                    "block_median_Ks",
                    "block_zscore",
                    "reason",
                ]
            )
            for event in events:
                writer.writerow(event.to_row())


def discover_species_directories(root: Path) -> Dict[str, Path]:
    """Discover MCScanX output directories under ``root``."""

    species_dirs: Dict[str, Path] = {}
    for path in sorted(Path(root).iterdir()):
        if not path.is_dir():
            continue
        gene_type_files = list(path.glob("*.gene_type"))
        if not gene_type_files:
            continue
        for file in gene_type_files:
            species = file.stem
            species_dirs[species] = path
    return species_dirs
