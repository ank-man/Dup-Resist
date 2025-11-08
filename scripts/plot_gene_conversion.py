#!/usr/bin/env python3
"""Produce visualizations summarizing gene conversion calls."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot gene conversion statistics.")
    parser.add_argument(
        "events",
        type=Path,
        nargs="+",
        help="Paths to *_gene_conversion.tsv files produced by Dup-Resist.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("gene_conversion_summary.png"),
        help="Path to the output figure (PNG).",
    )
    return parser.parse_args()


def load_events(paths: List[Path]) -> Dict[str, List[float]]:
    by_species: Dict[str, List[float]] = {}
    for path in paths:
        species = path.stem.replace("_gene_conversion", "")
        ks_values: List[float] = []
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                try:
                    ks = float(row["Ks"])
                except (KeyError, ValueError):
                    continue
                ks_values.append(ks)
        by_species[species] = ks_values
    return by_species


def main() -> None:
    args = parse_args()
    events = load_events(args.events)
    if not events:
        raise SystemExit("No gene conversion events were loaded.")

    fig, ax = plt.subplots(figsize=(8, 4 + len(events) * 0.5))
    species = list(events.keys())
    ks_series = [events[name] for name in species]
    ax.boxplot(ks_series, vert=False, labels=species, showfliers=True)
    ax.set_xlabel("Ks")
    ax.set_title("Distribution of Ks among putative gene conversion events")
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    print(f"Saved figure to {args.output}")


if __name__ == "__main__":
    main()
