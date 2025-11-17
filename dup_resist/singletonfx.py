"""SingletonFX pipeline for identifying conserved singleton functions.

This module implements the end-to-end workflow required by the ``singletonfx``
command line interface.  The pipeline extracts singleton genes from MCScanX
outputs, performs functional annotation with external tools (eggNOG-mapper
primary, InterProScan/Pfam optional), aggregates the annotations into a
long-format table, and carries out enrichment tests within and across genomes.

The implementation favours explicit, well documented helper functions so the
workflow can be reused programmatically and the external command invocations
are easy to audit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import math
import os
import sqlite3
import subprocess
import sys
import time
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

from matplotlib import pyplot as plt


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class TermKey:
    """Key identifying a functional term across genomes."""

    source: str
    term_type: str
    term_id: str
    term_name: str


@dataclass
class EnrichmentResult:
    """Fisher's exact test result for a term within a genome."""

    term: TermKey
    genome: str
    a_singleton_with: int
    b_singleton_without: int
    c_background_with: int
    d_background_without: int
    odds_ratio: float
    log2_odds_ratio: float
    p_value: float
    fdr: float
    significant: bool


AnnotationRecord = namedtuple(
    "AnnotationRecord",
    ["gene_id", "term_source", "term_type", "term_id", "term_name"],
)


def read_genome_list(path: Path) -> List[str]:
    """Return genome prefixes listed in ``path`` ignoring comments/empty lines."""

    genomes: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            clean = line.strip()
            if not clean or clean.startswith("#"):
                continue
            genomes.append(clean)
    if not genomes:
        raise ValueError(f"No genomes found in {path}")
    return genomes


def parse_gene_type_table(path: Path) -> Dict[str, int]:
    """Parse ``*.gene_type`` file mapping gene identifier to duplicate class."""

    mapping: Dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                LOGGER.debug("Skipping malformed gene_type line: %s", line)
                continue
            gene_id, class_str = fields[0], fields[1]
            try:
                mapping[gene_id] = int(class_str)
            except ValueError:
                LOGGER.debug("Non-integer duplicate class for %s: %s", gene_id, class_str)
    return mapping


def read_fasta(path: Path) -> Iterator[Tuple[str, str, str]]:
    """Yield (identifier, header, sequence) tuples from a FASTA file."""

    header: Optional[str] = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_chunks).replace("\n", "").strip()
                    identifier = header.split()[0][1:]
                    yield identifier, header.strip(), seq
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            seq = "".join(seq_chunks).replace("\n", "").strip()
            identifier = header.split()[0][1:]
            yield identifier, header.strip(), seq


def write_subset_fasta(source: Path, target: Path, keep_ids: Set[str]) -> int:
    """Write FASTA records whose identifiers are in ``keep_ids``; return count."""

    count = 0
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", encoding="utf-8") as in_handle, target.open(
        "w", encoding="utf-8"
    ) as out_handle:
        write_block = False
        for line in in_handle:
            if line.startswith(">"):
                identifier = line.split()[0][1:]
                write_block = identifier in keep_ids
                if write_block:
                    count += 1
                    out_handle.write(line)
            elif write_block:
                out_handle.write(line)
    return count


def count_fasta_records(path: Path) -> int:
    """Return the number of records in a FASTA file."""

    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def sha256sum(path: Path) -> str:
    """Compute the SHA-256 checksum of ``path``."""

    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_sqlite_cache(path: Path) -> sqlite3.Connection:
    """Initialise and return a connection to the metadata cache."""

    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS dataset_hash (
            genome TEXT NOT NULL,
            dataset TEXT NOT NULL,
            fasta_hash TEXT NOT NULL,
            n_records INTEGER NOT NULL,
            updated_at REAL NOT NULL,
            PRIMARY KEY (genome, dataset)
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS annotation_runs (
            genome TEXT NOT NULL,
            dataset TEXT NOT NULL,
            annotator TEXT NOT NULL,
            fasta_hash TEXT NOT NULL,
            status TEXT NOT NULL,
            updated_at REAL NOT NULL,
            PRIMARY KEY (genome, dataset, annotator)
        )
        """
    )
    conn.commit()
    return conn


def should_skip_annotation(
    conn: sqlite3.Connection,
    genome: str,
    dataset: str,
    annotator: str,
    fasta_hash: str,
    resume: bool,
) -> bool:
    """Return True when cached metadata allows skipping an annotation run."""

    if not resume:
        return False
    row = conn.execute(
        "SELECT fasta_hash, status FROM annotation_runs WHERE genome=? AND dataset=? AND annotator=?",
        (genome, dataset, annotator),
    ).fetchone()
    if row is None:
        return False
    cached_hash, status = row
    if status != "done":
        return False
    return cached_hash == fasta_hash


def record_annotation_status(
    conn: sqlite3.Connection,
    genome: str,
    dataset: str,
    annotator: str,
    fasta_hash: str,
    status: str,
) -> None:
    """Persist the completion status of an annotation run."""

    conn.execute(
        """
        INSERT INTO annotation_runs (genome, dataset, annotator, fasta_hash, status, updated_at)
        VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(genome, dataset, annotator)
        DO UPDATE SET fasta_hash=excluded.fasta_hash, status=excluded.status, updated_at=excluded.updated_at
        """,
        (genome, dataset, annotator, fasta_hash, status, time.time()),
    )
    conn.commit()


def update_dataset_hash(
    conn: sqlite3.Connection,
    genome: str,
    dataset: str,
    fasta_hash: str,
    n_records: int,
) -> None:
    """Store the checksum and record count for a FASTA dataset."""

    conn.execute(
        """
        INSERT INTO dataset_hash (genome, dataset, fasta_hash, n_records, updated_at)
        VALUES (?, ?, ?, ?, ?)
        ON CONFLICT(genome, dataset)
        DO UPDATE SET fasta_hash=excluded.fasta_hash, n_records=excluded.n_records, updated_at=excluded.updated_at
        """,
        (genome, dataset, fasta_hash, n_records, time.time()),
    )
    conn.commit()


def read_dataset_hash(
    conn: sqlite3.Connection, genome: str, dataset: str
) -> Optional[Tuple[str, int]]:
    """Return cached (hash, n_records) tuple for ``genome``/``dataset`` if present."""

    row = conn.execute(
        "SELECT fasta_hash, n_records FROM dataset_hash WHERE genome=? AND dataset=?",
        (genome, dataset),
    ).fetchone()
    if row is None:
        return None
    return str(row[0]), int(row[1])


def build_eggnog_command(
    fasta: Path, output_prefix: Path, threads: int, data_dir: Optional[Path]
) -> List[str]:
    cmd = [
        "emapper.py",
        "-i",
        str(fasta),
        "-o",
        str(output_prefix),
        "--cpu",
        str(threads),
    ]
    if data_dir is not None:
        cmd.extend(["--data_dir", str(data_dir)])
    cmd.append("--override")
    return cmd


def build_interpro_command(fasta: Path, output_tsv: Path, threads: int) -> List[str]:
    return [
        "interproscan.sh",
        "-i",
        str(fasta),
        "-f",
        "TSV",
        "-cpu",
        str(threads),
        "-o",
        str(output_tsv),
    ]


def build_pfam_command(
    fasta: Path, domtblout: Path, threads: int, database: Path
) -> List[str]:
    return [
        "hmmsearch",
        "--cpu",
        str(threads),
        "--domtblout",
        str(domtblout),
        str(database),
        str(fasta),
    ]


def run_command(cmd: Sequence[str], dry_run: bool) -> None:
    if dry_run:
        LOGGER.info("[dry-run] Would execute: %s", " ".join(cmd))
        return
    LOGGER.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_eggnog_annotations(path: Path) -> List[AnnotationRecord]:
    records: List[AnnotationRecord] = []
    if not path.exists():
        LOGGER.warning("eggNOG annotation file missing: %s", path)
        return records
    with path.open("r", encoding="utf-8") as handle:
        header_fields: Optional[List[str]] = None
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header_line = line.lstrip("#")
                if header_line.startswith("query") or header_line.startswith("query_name"):
                    header_fields = header_line.split("\t")
                continue
            if header_fields is None:
                LOGGER.warning("eggNOG header missing in %s", path)
                continue
            fields = line.split("\t")
            if len(fields) != len(header_fields):
                LOGGER.debug("Skipping malformed eggNOG line: %s", line)
                continue
            row = dict(zip(header_fields, fields))
            gene_id = row.get("query") or row.get("query_name") or row.get("#query")
            if not gene_id:
                continue
            description = row.get("Description", "")
            # GO terms
            go_terms = _split_terms(row.get("GOs"))
            for go in go_terms:
                records.append(
                    AnnotationRecord(gene_id, "eggNOG", "GO", go, "")
                )
            # KEGG Orthology
            for ko in _split_terms(row.get("KEGG_ko")):
                records.append(
                    AnnotationRecord(gene_id, "eggNOG", "KEGG_KO", ko, "")
                )
            for pathway in _split_terms(row.get("KEGG_Pathway")):
                records.append(
                    AnnotationRecord(gene_id, "eggNOG", "KEGG_Pathway", pathway, "")
                )
            for module in _split_terms(row.get("KEGG_Module")):
                records.append(
                    AnnotationRecord(gene_id, "eggNOG", "KEGG_Module", module, "")
                )
            for cog in _split_terms(row.get("COG_category")):
                records.append(
                    AnnotationRecord(gene_id, "eggNOG", "COG_Category", cog, "")
                )
            for og in _split_terms(row.get("eggNOG_OGs")):
                records.append(
                    AnnotationRecord(gene_id, "eggNOG", "eggNOG_OG", og, description)
                )
    return records


def parse_interproscan(path: Path) -> List[AnnotationRecord]:
    records: List[AnnotationRecord] = []
    if not path.exists():
        LOGGER.warning("InterProScan output missing: %s", path)
        return records
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields:
                continue
            if fields[0].startswith("#"):
                continue
            if len(fields) < 11:
                continue
            gene_id = fields[0]
            signature_acc = fields[4]
            signature_desc = fields[5]
            interpro_acc = fields[11] if len(fields) > 11 else ""
            interpro_desc = fields[12] if len(fields) > 12 else ""
            if signature_acc:
                records.append(
                    AnnotationRecord(
                        gene_id, "InterProScan", "Signature", signature_acc, signature_desc
                    )
                )
            if interpro_acc:
                records.append(
                    AnnotationRecord(
                        gene_id, "InterProScan", "InterPro", interpro_acc, interpro_desc
                    )
                )
            if len(fields) > 13 and fields[13]:
                for go in _split_terms(fields[13]):
                    records.append(
                        AnnotationRecord(gene_id, "InterProScan", "GO", go, "")
                    )
            if len(fields) > 14 and fields[14]:
                for pathway in _split_terms(fields[14]):
                    records.append(
                        AnnotationRecord(
                            gene_id, "InterProScan", "Pathway", pathway, ""
                        )
                    )
    return records


def parse_pfam_domtbl(path: Path) -> List[AnnotationRecord]:
    records: List[AnnotationRecord] = []
    if not path.exists():
        LOGGER.warning("Pfam domtblout missing: %s", path)
        return records
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            target_name = fields[0]
            hmm_name = fields[3]
            hmm_acc = fields[4]
            term_id = hmm_acc.split(".")[0] if hmm_acc else hmm_name
            term_name = hmm_name
            records.append(
                AnnotationRecord(target_name, "Pfam", "Pfam", term_id, term_name)
            )
    return records


def _split_terms(value: Optional[str]) -> List[str]:
    if not value or value in {"-", ""}:
        return []
    separators = [",", "|", ";"]
    for sep in separators:
        if sep in value:
            parts = [part.strip() for part in value.split(sep)]
            return [part for part in parts if part]
    return [value.strip()]


class GoSlimMapper:
    """Map GO terms to goslim_generic categories."""

    def __init__(self, obo_path: Path) -> None:
        self.obo_path = obo_path
        self.parents: Dict[str, List[str]] = defaultdict(list)
        self.slim_terms: Set[str] = set()
        self._load()
        self._memo: Dict[str, Set[str]] = {}

    def _load(self) -> None:
        current_id: Optional[str] = None
        with self.obo_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line == "[Term]":
                    current_id = None
                    continue
                if not line:
                    continue
                if line.startswith("id: "):
                    current_id = line.split("id: ", 1)[1]
                elif line.startswith("is_a: ") and current_id:
                    parent = line.split("is_a: ", 1)[1].split()[0]
                    self.parents[current_id].append(parent)
                elif line.startswith("subset: ") and current_id:
                    subset = line.split("subset: ", 1)[1]
                    if subset == "goslim_generic":
                        self.slim_terms.add(current_id)

    def map(self, term: str) -> Set[str]:
        if term in self._memo:
            return self._memo[term]
        if term in self.slim_terms:
            self._memo[term] = {term}
            return {term}
        seen: Set[str] = set()
        queue: List[str] = [term]
        results: Set[str] = set()
        while queue:
            current = queue.pop()
            if current in seen:
                continue
            seen.add(current)
            if current in self.slim_terms:
                results.add(current)
                continue
            for parent in self.parents.get(current, []):
                queue.append(parent)
        if not results:
            results = {term}
        self._memo[term] = results
        return results


def benjamini_hochberg(p_values: List[float]) -> List[float]:
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, (idx, pval) in enumerate(indexed, start=1):
        adj = min(prev, pval * n / rank)
        adjusted[idx] = adj
        prev = adj
    # ensure monotonicity by iterating backwards
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    return adjusted


def fisher_exact(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    """Compute odds ratio and two-sided Fisher exact p-value for a 2x2 table."""

    odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    n = a + b + c + d
    row_sum = a + b
    col_sum = a + c
    min_a = max(0, row_sum - (n - col_sum))
    max_a = min(row_sum, col_sum)

    def hypergeom(x: int) -> float:
        return (
            math.comb(row_sum, x)
            * math.comb(n - row_sum, col_sum - x)
            / math.comb(n, col_sum)
        )

    obs_prob = hypergeom(a)
    p_value = 0.0
    for x in range(min_a, max_a + 1):
        prob = hypergeom(x)
        if prob <= obs_prob + 1e-12:
            p_value += prob
    p_value = min(p_value, 1.0)
    return odds_ratio, p_value


def binomial_test_two_sided(k: int, n: int, p: float = 0.5) -> float:
    if n == 0:
        return 1.0
    prob_k = math.comb(n, k) * (p**k) * ((1 - p) ** (n - k))
    total = 0.0
    for i in range(0, n + 1):
        prob = math.comb(n, i) * (p**i) * ((1 - p) ** (n - i))
        if prob <= prob_k + 1e-12:
            total += prob
    return min(1.0, total)


class SingletonFXPipeline:
    """Coordinator for the singleton functional enrichment workflow."""

    def __init__(
        self,
        genomes: Sequence[str],
        pep_dir: Path,
        mcscanx_dir: Path,
        out_dir: Path,
        annotators: Sequence[str],
        threads: int,
        min_genes_per_term: int,
        alpha: float,
        resume: bool,
        dry_run: bool,
        go_slim_mapper: Optional[GoSlimMapper] = None,
        eggnog_data_dir: Optional[Path] = None,
        pfam_db: Optional[Path] = None,
        plot_top_n: int = 20,
    ) -> None:
        self.genomes = list(genomes)
        self.pep_dir = pep_dir
        self.mcscanx_dir = mcscanx_dir
        self.out_dir = out_dir
        self.annotators = [annotator.lower() for annotator in annotators]
        if not self.annotators:
            raise ValueError("At least one annotator must be specified")
        self.threads = threads
        self.min_genes_per_term = min_genes_per_term
        self.alpha = alpha
        self.resume = resume
        self.dry_run = dry_run
        self.go_slim_mapper = go_slim_mapper
        self.eggnog_data_dir = eggnog_data_dir
        self.pfam_db = pfam_db
        self.plot_top_n = max(1, plot_top_n)

        self.cache_path = self.out_dir / "singletonfx_cache.sqlite"
        self.conn = ensure_sqlite_cache(self.cache_path)
        self.per_genome_enrichment: Dict[str, List[EnrichmentResult]] = {}

    def run(self) -> None:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        for genome in self.genomes:
            try:
                self._process_genome(genome)
            except Exception as exc:  # pragma: no cover - keep pipeline resilient
                LOGGER.exception("Failed to process %s: %s", genome, exc)
        self._write_cross_species_outputs()
        self.conn.close()

    # ------------------------------------------------------------------
    # Genome level processing
    # ------------------------------------------------------------------
    def _process_genome(self, genome: str) -> None:
        LOGGER.info("Processing genome %s", genome)
        pep_path = self.pep_dir / f"{genome}.pep.fa"
        gene_type_path = self._resolve_gene_type(genome)
        out_genome_dir = self.out_dir / genome
        out_genome_dir.mkdir(parents=True, exist_ok=True)
        if not pep_path.exists():
            raise FileNotFoundError(f"Peptide FASTA missing for {genome}: {pep_path}")
        if not gene_type_path.exists():
            raise FileNotFoundError(f"gene_type file missing for {genome}: {gene_type_path}")

        gene_classes = parse_gene_type_table(gene_type_path)
        singleton_ids = {gid for gid, cls in gene_classes.items() if cls == 0}
        if not singleton_ids:
            LOGGER.warning("No singleton genes found for %s; skipping", genome)
            return

        singleton_fasta = out_genome_dir / f"{genome}.singletons.fa"
        n_singleton_records = self._write_singleton_fasta(
            genome, pep_path, singleton_fasta, singleton_ids
        )

        dataset_hashes: Dict[str, Tuple[str, int]] = {}
        if self.dry_run:
            cached_hash = read_dataset_hash(self.conn, genome, "singletons")
            singleton_hash = cached_hash[0] if cached_hash else "dryrun"
        else:
            singleton_hash = sha256sum(singleton_fasta)
            update_dataset_hash(
                self.conn, genome, "singletons", singleton_hash, n_singleton_records
            )
        dataset_hashes["singletons"] = (singleton_hash, n_singleton_records)

        total_records = count_fasta_records(pep_path)
        all_hash = sha256sum(pep_path)
        dataset_hashes["all"] = (all_hash, total_records)
        if not self.dry_run:
            update_dataset_hash(self.conn, genome, "all", all_hash, total_records)

        annotations_dir = out_genome_dir / "annotations"
        annotations_dir.mkdir(exist_ok=True)

        for dataset, (fasta_hash, _) in dataset_hashes.items():
            fasta_path = singleton_fasta if dataset == "singletons" else pep_path
            for annotator in self.annotators:
                self._run_annotator(
                    genome,
                    dataset,
                    annotator,
                    fasta_path,
                    fasta_hash,
                    annotations_dir,
                )

        singleton_annotations = self._load_annotations(genome, "singletons", annotations_dir)
        all_annotations = self._load_annotations(genome, "all", annotations_dir)

        func_path = out_genome_dir / f"{genome}.func.tsv"
        self._write_annotation_table(func_path, singleton_annotations)

        enrichment = self._compute_enrichment(
            genome,
            singleton_ids,
            singleton_annotations,
            gene_classes,
            all_annotations,
        )
        self.per_genome_enrichment[genome] = enrichment
        enrich_path = out_genome_dir / f"{genome}.singleton_enrichment.tsv"
        self._write_enrichment_table(enrich_path, enrichment)
        plot_path = out_genome_dir / f"{genome}.singleton_enrichment.png"
        self._plot_enrichment(plot_path, genome, enrichment)

    def _write_singleton_fasta(
        self, genome: str, pep_path: Path, singleton_fasta: Path, singleton_ids: Set[str]
    ) -> int:
        cached = read_dataset_hash(self.conn, genome, "singletons")
        if self.resume and singleton_fasta.exists():
            current_hash = sha256sum(singleton_fasta)
            if cached and cached[0] == current_hash:
                LOGGER.info("Reusing existing singleton FASTA for %s", pep_path.stem)
                return cached[1]
        if self.dry_run:
            LOGGER.info("[dry-run] Would write singleton FASTA: %s", singleton_fasta)
            n_records = len(singleton_ids)
        else:
            n_records = write_subset_fasta(pep_path, singleton_fasta, singleton_ids)
            LOGGER.info("Wrote %d singleton peptide records to %s", n_records, singleton_fasta)
        return n_records

    def _resolve_gene_type(self, genome: str) -> Path:
        direct = self.mcscanx_dir / genome / f"{genome}.gene_type"
        if direct.exists():
            return direct
        fallback = self.mcscanx_dir / f"{genome}.gene_type"
        return fallback

    def _run_annotator(
        self,
        genome: str,
        dataset: str,
        annotator: str,
        fasta_path: Path,
        fasta_hash: str,
        annotations_dir: Path,
    ) -> None:
        annotator = annotator.lower()
        annot_dir = annotations_dir / annotator
        annot_dir.mkdir(parents=True, exist_ok=True)

        if should_skip_annotation(
            self.conn, genome, dataset, annotator, fasta_hash, self.resume
        ):
            LOGGER.info(
                "Skipping %s annotation for %s (%s) due to cached hash",
                annotator,
                genome,
                dataset,
            )
            return

        status = "done"
        try:
            if annotator == "eggnog":
                output_prefix = annot_dir / f"{genome}_{dataset}"
                cmd = build_eggnog_command(
                    fasta_path,
                    output_prefix,
                    self.threads,
                    self.eggnog_data_dir,
                )
                run_command(cmd, self.dry_run)
            elif annotator == "interproscan":
                output_tsv = annot_dir / f"{genome}_{dataset}.interproscan.tsv"
                cmd = build_interpro_command(fasta_path, output_tsv, self.threads)
                run_command(cmd, self.dry_run)
            elif annotator == "pfam":
                if self.pfam_db is None:
                    raise ValueError("Pfam database path must be provided via --pfam-db")
                domtblout = annot_dir / f"{genome}_{dataset}.pfam.domtblout"
                cmd = build_pfam_command(
                    fasta_path, domtblout, self.threads, self.pfam_db
                )
                run_command(cmd, self.dry_run)
            else:
                raise ValueError(f"Unsupported annotator: {annotator}")
        except Exception:
            status = "failed"
            record_annotation_status(
                self.conn, genome, dataset, annotator, fasta_hash, status
            )
            raise
        else:
            final_status = "dry_run" if self.dry_run else status
            record_annotation_status(
                self.conn, genome, dataset, annotator, fasta_hash, final_status
            )

    def _load_annotations(
        self, genome: str, dataset: str, annotations_dir: Path
    ) -> List[AnnotationRecord]:
        records: List[AnnotationRecord] = []
        for annotator in self.annotators:
            annot_dir = annotations_dir / annotator
            if annotator == "eggnog":
                path = annot_dir / f"{genome}_{dataset}.emapper.annotations"
                records.extend(parse_eggnog_annotations(path))
            elif annotator == "interproscan":
                path = annot_dir / f"{genome}_{dataset}.interproscan.tsv"
                records.extend(parse_interproscan(path))
            elif annotator == "pfam":
                path = annot_dir / f"{genome}_{dataset}.pfam.domtblout"
                records.extend(parse_pfam_domtbl(path))
        if self.go_slim_mapper:
            expanded: List[AnnotationRecord] = []
            for rec in records:
                if rec.term_type != "GO":
                    expanded.append(rec)
                    continue
                for mapped in self.go_slim_mapper.map(rec.term_id):
                    expanded.append(
                        AnnotationRecord(
                            rec.gene_id,
                            rec.term_source,
                            "GO_slim",
                            mapped,
                            rec.term_name,
                        )
                    )
            records = expanded
        return records

    def _write_annotation_table(
        self, path: Path, annotations: List[AnnotationRecord]
    ) -> None:
        if self.dry_run:
            LOGGER.info("[dry-run] Would write annotation table: %s", path)
            return
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                ["gene_id", "term_source", "term_type", "term_id", "term_name"]
            )
            for rec in annotations:
                writer.writerow(rec)

    def _compute_enrichment(
        self,
        genome: str,
        singleton_ids: Set[str],
        singleton_annotations: List[AnnotationRecord],
        gene_classes: Dict[str, int],
        all_annotations: List[AnnotationRecord],
    ) -> List[EnrichmentResult]:
        singleton_terms: Dict[TermKey, Set[str]] = defaultdict(set)
        for rec in singleton_annotations:
            term = TermKey(rec.term_source, rec.term_type, rec.term_id, rec.term_name)
            singleton_terms[term].add(rec.gene_id)

        background_terms: Dict[TermKey, Set[str]] = defaultdict(set)
        for rec in all_annotations:
            term = TermKey(rec.term_source, rec.term_type, rec.term_id, rec.term_name)
            background_terms[term].add(rec.gene_id)

        singleton_total = len(singleton_ids)
        background_ids = {gid for gid, cls in gene_classes.items() if cls != 0}
        background_total = len(background_ids)

        results: List[EnrichmentResult] = []
        all_terms = set(background_terms.keys()) | set(singleton_terms.keys())
        for term in sorted(all_terms, key=lambda t: (t.source, t.term_type, t.term_id)):
            singleton_with = len(singleton_terms.get(term, set()))
            if singleton_with < self.min_genes_per_term:
                continue
            background_with = len(
                {gid for gid in background_terms.get(term, set()) if gid in background_ids}
            )
            a = singleton_with
            b = singleton_total - singleton_with
            c = background_with
            d = background_total - background_with
            if (a + b + c + d) == 0:
                continue
            odds_ratio, p_value = fisher_exact(a, b, c, d)
            log2_or = math.log2(odds_ratio) if odds_ratio > 0 else float("nan")
            results.append(
                EnrichmentResult(
                    term=term,
                    genome=genome,
                    a_singleton_with=a,
                    b_singleton_without=b,
                    c_background_with=c,
                    d_background_without=d,
                    odds_ratio=odds_ratio,
                    log2_odds_ratio=log2_or,
                    p_value=p_value,
                    fdr=1.0,
                    significant=False,
                )
            )

        p_values = [res.p_value for res in results]
        fdrs = benjamini_hochberg(p_values)
        for res, fdr in zip(results, fdrs):
            res.fdr = fdr
            res.significant = fdr <= self.alpha
        return results

    def _write_enrichment_table(
        self, path: Path, results: List[EnrichmentResult]
    ) -> None:
        if self.dry_run:
            LOGGER.info("[dry-run] Would write enrichment table: %s", path)
            return
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "term_source",
                    "term_type",
                    "term_id",
                    "term_name",
                    "singleton_with_term",
                    "singleton_without_term",
                    "background_with_term",
                    "background_without_term",
                    "odds_ratio",
                    "log2_odds_ratio",
                    "p_value",
                    "fdr",
                    "significant",
                ]
            )
            for res in results:
                writer.writerow(
                    [
                        res.term.source,
                        res.term.term_type,
                        res.term.term_id,
                        res.term.term_name,
                        res.a_singleton_with,
                        res.b_singleton_without,
                        res.c_background_with,
                        res.d_background_without,
                        f"{res.odds_ratio:.6g}",
                        f"{res.log2_odds_ratio:.6g}",
                        f"{res.p_value:.6g}",
                        f"{res.fdr:.6g}",
                        "yes" if res.significant else "no",
                    ]
                )

    def _plot_enrichment(
        self, path: Path, genome: str, results: List[EnrichmentResult]
    ) -> None:
        if self.dry_run:
            LOGGER.info("[dry-run] Would write enrichment plot: %s", path)
            return
        filtered = [res for res in results if math.isfinite(res.log2_odds_ratio)]
        if not filtered:
            LOGGER.info("No enrichment results to plot for %s", genome)
            return

        sorted_results = sorted(
            filtered,
            key=lambda res: (res.fdr, -res.log2_odds_ratio if math.isfinite(res.log2_odds_ratio) else 0.0),
        )
        top_results = sorted_results[: self.plot_top_n]
        if not top_results:
            LOGGER.info("No enrichment results remain after filtering for %s", genome)
            return

        labels: List[str] = []
        values: List[float] = []
        colors: List[str] = []
        from matplotlib.patches import Patch

        for res in top_results:
            term_label = f"{res.term.source}:{res.term.term_id}"
            if res.term.term_name and res.term.term_name.strip():
                term_label = f"{term_label} — {res.term.term_name.strip()}"
            labels.append(term_label)
            values.append(res.log2_odds_ratio)
            colors.append("#2166ac" if res.significant else "#bdbdbd")

        height = max(3.5, 0.4 * len(values) + 1.5)
        fig, ax = plt.subplots(figsize=(10, height))
        y_positions = range(len(values))
        ax.barh(y_positions, values, color=colors, edgecolor="#333333", linewidth=0.6)
        ax.axvline(0.0, color="#636363", linestyle="--", linewidth=1.0)
        ax.set_yticks(list(y_positions))
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel("log2 odds ratio (Singleton vs Background)")
        ax.set_title(
            f"{genome} singleton functional enrichment (top {len(values)} terms)"
        )

        legend_handles = [
            Patch(facecolor="#2166ac", edgecolor="#333333", label=f"FDR ≤ {self.alpha}"),
            Patch(facecolor="#bdbdbd", edgecolor="#333333", label="Not significant"),
        ]
        ax.legend(handles=legend_handles, loc="lower right", frameon=False)
        ax.grid(axis="x", color="#e0e0e0", linestyle="-", linewidth=0.5, alpha=0.6)
        fig.tight_layout()
        try:
            fig.savefig(path, dpi=300)
            LOGGER.info("Wrote enrichment plot to %s", path)
        except Exception:
            LOGGER.exception("Failed to save enrichment plot for %s", genome)
        finally:
            plt.close(fig)

    # ------------------------------------------------------------------
    # Cross genome summaries
    # ------------------------------------------------------------------
    def _write_cross_species_outputs(self) -> None:
        if not self.per_genome_enrichment:
            LOGGER.warning("No enrichment results available for cross-species summary")
            return

        summary_path = self.out_dir / "CONSISTENT_SINGLETON_CLASSES.tsv"
        matrix_path = self.out_dir / "SINGLETON_LOG2OR_MATRIX.tsv"

        cross_stats = self._compute_cross_species_stats()
        if self.dry_run:
            LOGGER.info("[dry-run] Would write cross-species summary: %s", summary_path)
            LOGGER.info("[dry-run] Would write log2 OR matrix: %s", matrix_path)
            return

        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "term_source",
                    "term_type",
                    "term_id",
                    "term_name",
                    "n_genomes_tested",
                    "n_significant",
                    "fraction_significant",
                    "median_odds_ratio",
                    "mantel_haenszel_or",
                    "sign_test_p",
                    "q_sign",
                ]
            )
            for stats in cross_stats:
                writer.writerow(
                    [
                        stats["term"].source,
                        stats["term"].term_type,
                        stats["term"].term_id,
                        stats["term"].term_name,
                        stats["n_tested"],
                        stats["n_significant"],
                        f"{stats['fraction_significant']:.6g}",
                        f"{stats['median_or']:.6g}",
                        f"{stats['mantel_haenszel_or']:.6g}",
                        f"{stats['sign_test_p']:.6g}",
                        f"{stats['q_sign']:.6g}",
                    ]
                )

        genomes = sorted(self.per_genome_enrichment.keys())
        terms = [stats["term"] for stats in cross_stats]
        with matrix_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            header = ["term_source", "term_type", "term_id", "term_name"] + genomes
            writer.writerow(header)
            for term in terms:
                row = [term.source, term.term_type, term.term_id, term.term_name]
                for genome in genomes:
                    log2_or = self._lookup_log2_or(genome, term)
                    row.append(f"{log2_or:.6g}" if not math.isnan(log2_or) else "NA")
                writer.writerow(row)

    def _lookup_log2_or(self, genome: str, term: TermKey) -> float:
        results = self.per_genome_enrichment.get(genome, [])
        for res in results:
            if res.term == term:
                return res.log2_odds_ratio
        return float("nan")

    def _compute_cross_species_stats(self) -> List[Dict[str, object]]:
        term_to_tables: Dict[TermKey, List[Tuple[int, int, int, int, EnrichmentResult]]] = defaultdict(list)
        for genome, results in self.per_genome_enrichment.items():
            for res in results:
                term_to_tables[res.term].append(
                    (
                        res.a_singleton_with,
                        res.b_singleton_without,
                        res.c_background_with,
                        res.d_background_without,
                        res,
                    )
                )

        cross_stats: List[Dict[str, object]] = []
        sign_p_values: List[float] = []
        term_order: List[TermKey] = []
        for term in sorted(term_to_tables.keys(), key=lambda t: (t.source, t.term_type, t.term_id)):
            tables = term_to_tables[term]
            n_tested = len(tables)
            n_significant = sum(1 for _, _, _, _, res in tables if res.significant)
            frac_significant = n_significant / n_tested if n_tested else 0.0
            odds_ratios = [res.odds_ratio for _, _, _, _, res in tables if math.isfinite(res.odds_ratio)]
            med_or = median(odds_ratios) if odds_ratios else float("nan")
            mh_or = self._mantel_haenszel_or(tables)
            n_positive = sum(1 for _, _, _, _, res in tables if res.odds_ratio > 1.0)
            sign_p = binomial_test_two_sided(n_positive, n_tested, 0.5)
            sign_p_values.append(sign_p)
            term_order.append(term)
            cross_stats.append(
                {
                    "term": term,
                    "n_tested": n_tested,
                    "n_significant": n_significant,
                    "fraction_significant": frac_significant,
                    "median_or": med_or,
                    "mantel_haenszel_or": mh_or,
                    "sign_test_p": sign_p,
                    "q_sign": 1.0,
                }
            )

        q_values = benjamini_hochberg(sign_p_values)
        for stats, q in zip(cross_stats, q_values):
            stats["q_sign"] = q
        return cross_stats

    def _mantel_haenszel_or(
        self, tables: List[Tuple[int, int, int, int, EnrichmentResult]]
    ) -> float:
        num = 0.0
        den = 0.0
        for a, b, c, d, _ in tables:
            n = a + b + c + d
            if n == 0:
                continue
            num += (a * d) / n
            den += (b * c) / n
        if den == 0.0:
            return float("nan")
        return num / den


def build_singletonfx_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Identify conserved singleton gene functions across genomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--genome-list",
        type=Path,
        required=True,
        help="File containing genome prefixes (one per line)",
    )
    parser.add_argument(
        "--pep-dir",
        type=Path,
        default=Path("/home/ankush/PGDD/files/gold/mcscanx_out/pep/processed"),
        help="Directory containing GEN.pep.fa peptide FASTAs",
    )
    parser.add_argument(
        "--mcscanx-dir",
        type=Path,
        required=True,
        help="Directory with MCScanX outputs (GEN/GEN.gene_type)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("singletonfx_out"),
        help="Output directory for all derived files",
    )
    parser.add_argument(
        "--annotators",
        type=str,
        default="eggnog",
        help="Comma separated annotators to run (eggnog, interproscan, pfam)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Worker threads for external tools",
    )
    parser.add_argument(
        "--min-genes-per-term",
        type=int,
        default=5,
        help="Minimum singleton genes annotated with a term to test enrichment",
    )
    parser.add_argument(
        "--plot-top-n",
        type=int,
        default=20,
        help="Number of top enriched terms to visualise per genome",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR threshold for significance",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse cached results when FASTA hashes match",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without running external tools or writing files",
    )
    parser.add_argument(
        "--go-slim-obo",
        type=Path,
        help="Optional path to goslim_generic.obo for GO term collapsing",
    )
    parser.add_argument(
        "--eggnog-data-dir",
        type=Path,
        default=Path.home() / "eggnog_data",
        help="eggNOG-mapper data directory",
    )
    parser.add_argument(
        "--pfam-db",
        type=Path,
        help="Pfam-A HMM database file for hmmsearch",
    )
    return parser


def run_singletonfx(args: argparse.Namespace) -> None:
    annotators = [part.strip() for part in args.annotators.split(",") if part.strip()]
    genomes = read_genome_list(args.genome_list)
    go_mapper = GoSlimMapper(args.go_slim_obo) if args.go_slim_obo else None

    pipeline = SingletonFXPipeline(
        genomes=genomes,
        pep_dir=args.pep_dir,
        mcscanx_dir=args.mcscanx_dir,
        out_dir=args.out_dir,
        annotators=annotators,
        threads=args.threads,
        min_genes_per_term=args.min_genes_per_term,
        alpha=args.alpha,
        resume=args.resume,
        dry_run=args.dry_run,
        go_slim_mapper=go_mapper,
        eggnog_data_dir=args.eggnog_data_dir,
        pfam_db=args.pfam_db,
        plot_top_n=args.plot_top_n,
    )
    pipeline.run()


__all__ = [
    "SingletonFXPipeline",
    "build_singletonfx_parser",
    "run_singletonfx",
]

