"""Command line entry point for the singletonfx workflow."""

from __future__ import annotations

import logging

from dup_resist.singletonfx import build_singletonfx_parser, run_singletonfx


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s - %(message)s",
    )


def main() -> None:
    configure_logging()
    parser = build_singletonfx_parser()
    args = parser.parse_args()
    run_singletonfx(args)


if __name__ == "__main__":  # pragma: no cover
    main()

