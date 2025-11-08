#!/usr/bin/env python3
"""CLI entry point for running Dup-Resist gene conversion detection."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

from dup_resist.gene_conversion import GeneConversionDetector, discover_species_directories


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Detect putative gene conversion events among WGD duplicates."
    )
    parser.add_argument(
        "mcscanx_root",
        type=Path,
        help="Path to the directory containing per-species MCScanX outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("gene_conversion_results"),
        help="Directory where gene conversion tables will be written.",
    )
    parser.add_argument(
        "--duplicate-class",
        default="WGD",
        help="Duplicate class to analyze (default: WGD).",
    )
    parser.add_argument(
        "--min-block-size",
        type=int,
        default=4,
        help="Minimum number of duplicate pairs in a block to analyze.",
    )
    parser.add_argument(
        "--min-pairs-with-ks",
        type=int,
        default=3,
        help="Minimum number of WGD pairs with Ks estimates to consider a block.",
    )
    parser.add_argument(
        "--z-threshold",
        type=float,
        default=1.5,
        help="Z-score threshold for marking unusually small Ks values.",
    )
    parser.add_argument(
        "--absolute-ks-threshold",
        type=float,
        default=None,
        help="Optional absolute Ks cutoff for calling gene conversion.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    detector = GeneConversionDetector(
        duplicate_class=args.duplicate_class,
        min_block_size=args.min_block_size,
        z_threshold=args.z_threshold,
        absolute_ks_threshold=args.absolute_ks_threshold,
        min_pairs_with_ks=args.min_pairs_with_ks,
    )

    species_dirs = discover_species_directories(args.mcscanx_root)
    if not species_dirs:
        raise SystemExit("No MCScanX species directories found.")

    summary_rows: List[Tuple[str, int]] = []
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "summary.tsv"

    for species, directory in species_dirs.items():
        events = detector.from_species_directory(species, directory, args.output_dir)
        summary_rows.append((species, len(events)))
        print(f"Processed {species}: {len(events)} events")

    with summary_path.open("w") as handle:
        handle.write("species\tevent_count\n")
        for species, count in summary_rows:
            handle.write(f"{species}\t{count}\n")


if __name__ == "__main__":
    main()
