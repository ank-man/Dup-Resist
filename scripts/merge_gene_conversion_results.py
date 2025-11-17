#!/usr/bin/env python3
"""Merge per-species gene conversion calls into a single table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List, Optional


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge Dup-Resist event tables.")
    parser.add_argument(
        "inputs",
        type=Path,
        nargs="+",
        help="Paths to *_gene_conversion.tsv files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("all_gene_conversion_events.tsv"),
        help="Output table path.",
    )
    return parser.parse_args()


def merge_tables(inputs: Iterable[Path], output: Path) -> None:
    header: Optional[List[str]] = None
    with output.open("w", newline="") as out_handle:
        writer = None
        for path in inputs:
            with path.open() as handle:
                reader = csv.reader(handle, delimiter="\t")
                rows = list(reader)
                if not rows:
                    continue
                if header is None:
                    header = rows[0]
                    writer = csv.writer(out_handle, delimiter="\t")
                    writer.writerow(header)
                for row in rows[1:]:
                    writer.writerow(row)
    if header is None:
        raise SystemExit("No rows were merged; please verify the input paths.")


def main() -> None:
    args = parse_args()
    merge_tables(args.inputs, args.output)
    print(f"Merged {len(args.inputs)} tables into {args.output}")


if __name__ == "__main__":
    main()
