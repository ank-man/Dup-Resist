# Dup-Resist

Dup-Resist is a small collection of Python utilities for analysing duplicated
gene repertoires.  It contains workflows for detecting putative gene conversion
events among whole genome duplication (WGD) pairs identified by MCScanX and for
profiling functionally conserved singleton (duplicate-resistant) genes across
multiple genomes.

## Functional singleton discovery with `singletonfx`

`singletonfx` automates the identification and functional annotation of
singleton genes (duplicate class 0) across many MCScanX runs.  It extracts the
singleton peptide sequences, annotates them via eggNOG-mapper (with optional
InterProScan/Pfam fallbacks), performs per-genome Fisher's exact tests against
the full gene set, and aggregates consistent enrichments across species.

### Expected inputs

- A text file listing genome prefixes (one per line), e.g. `genomes.txt`.
- MCScanX output directories containing `GEN/GEN.gene_type` and related files.
- Peptide FASTAs named `GEN.pep.fa` located under a peptide directory (default:
  `/home/ankush/PGDD/files/gold/mcscanx_out/pep/processed`).

### Quick start

```bash
singletonfx \
  --genome-list genomes.txt \
  --mcscanx-dir /path/to/mcscanx_intra/species \
  --pep-dir /home/ankush/PGDD/files/gold/mcscanx_out/pep/processed \
  --out-dir singletonfx_results \
  --threads 16 \
  --annotators eggnog,interproscan
```

The pipeline is incremental by default: it stores SHA-256 hashes of the FASTA
subsets in `singletonfx_cache.sqlite` and skips annotation runs whose inputs have
not changed when invoked with `--resume`.  Use `--dry-run` to inspect actions
without executing external tools.  Optional switches include `--go-slim-obo` to
collapse GO terms to the goslim_generic subset, `--pfam-db` to point at a Pfam-A
HMM database for the Pfam annotator, and `--plot-top-n` to control how many
terms appear in the per-genome enrichment figures.

### Outputs

For each genome `GEN` the command writes the following under
`OUT_DIR/GEN/`:

- `GEN.singletons.fa` – singleton peptide subset used for annotation.
- `GEN.func.tsv` – long-format table with columns
  `gene_id, term_source, term_type, term_id, term_name` merging all annotation
  providers.
- `GEN.singleton_enrichment.tsv` – Fisher's exact test results (odds ratios,
  p-values, FDRs, significance flags) comparing singleton genes against the
  remaining gene set.
- `GEN.singleton_enrichment.png` – horizontal bar plot showing the top enriched
  terms coloured by significance (FDR ≤ α vs. not significant).
- `annotations/` – raw outputs from eggNOG-mapper, InterProScan, and/or Pfam.

Once all genomes have been processed two cross-species summaries appear in the
output directory:

- `CONSISTENT_SINGLETON_CLASSES.tsv` – statistics describing how often a term is
  enriched across genomes, including Mantel–Haenszel common odds ratios and
  binomial sign tests (with FDR correction).
- `SINGLETON_LOG2OR_MATRIX.tsv` – heatmap-ready matrix of per-genome log2 odds
  ratios.

## Gene conversion detection utilities

The gene conversion detector expects that MCScanX has already been executed with
`duplicate_gene_classifier` and `add_ka_and_ks_to_collinearity.pl`, producing the
following files for each species:

- `*.collinearity`
- `*.gene_type`
- `*.kaks.tsv`

The scripts included here reimplement the functionality of the gene conversion
workflow from [Scripts_for_GB](https://github.com/qiao-xin/Scripts_for_GB/tree/master/detect_gene_conversion)
while adapting it to a multi-species directory layout.

## Repository layout

```
mcscanx_intra/
    species_a/
        species_a.collinearity
        species_a.gene_type
        species_a.kaks.tsv
    species_b/
        ...
gene_conversion_results/
```

All scripts reside in the `scripts/` directory and rely on the reusable
components that live inside the `dup_resist/` package.

## Installation

The utilities only depend on the Python standard library with the exception of
`matplotlib` for plotting. Install it with pip if necessary:

```bash
python -m pip install matplotlib
```

For convenience you may place the repository on your `PYTHONPATH` or install it
in editable mode:

```bash
python -m pip install -e .
```

## Detecting gene conversion events

To scan every species directory under `mcscanx_intra/` and collect candidate gene
conversion events among WGD pairs:

```bash
python scripts/detect_gene_conversion.py mcscanx_intra --output-dir gene_conversion_results
```

Key options:

- `--duplicate-class`: choose a different duplicate class (default: `WGD`).
- `--min-block-size`: minimum number of duplicate pairs per collinear block to
  evaluate (default: 4).
- `--min-pairs-with-ks`: minimum number of WGD pairs in a block that must have a
  Ks estimate (default: 3).
- `--z-threshold`: z-score cutoff for labeling an unusually small Ks (default:
  1.5).
- `--absolute-ks-threshold`: optional absolute Ks threshold; any pair with Ks
  below this value is marked as a gene conversion candidate.

The script writes one `*_gene_conversion.tsv` file per species alongside an
aggregated `summary.tsv` listing the number of events per species.

## Merging results across species

Once individual species have been processed you can combine the event tables
into a single summary:

```bash
python scripts/merge_gene_conversion_results.py gene_conversion_results/*_gene_conversion.tsv --output all_gene_conversion_events.tsv
```

## Plotting

To visualize the Ks distribution of the detected events across species, use the
plotting helper:

```bash
python scripts/plot_gene_conversion.py gene_conversion_results/*_gene_conversion.tsv --output gene_conversion_summary.png
```

The resulting PNG contains a horizontal boxplot for each species, providing a
quick overview of the Ks ranges among candidate conversion events.

## Reproducibility notes

- All thresholds are configurable through CLI options.
- The detector only considers pairs where both genes belong to the chosen
  duplicate class (WGD by default).
- Blocks must have a sufficient number of Ks estimates before they are
  evaluated.
- Candidate events are flagged when their Ks z-score is more than the specified
  threshold below the block median, optionally combined with an absolute Ks
  cutoff.

## Extending the toolkit

The `dup_resist` package exposes a programmatic API that can be reused in custom
workflows. See `dup_resist/gene_conversion.py` for details on the
`GeneConversionDetector` class and helper utilities.
