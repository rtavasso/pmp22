# Project A — PMP22 regulatory atlas

This repository is an executable, audit-first implementation of the research
contract in `project_a_regulatory_atlas.html`.  It deliberately separates
measured evidence, model predictions, and proposed experiments.  The included
records are a curated **starting atlas**, not new biological or clinical
findings.

**Read the succinct findings report:** [`report.html`](report.html)

## Quick start

```bash
python -m pmp22_atlas validate
python -m pmp22_atlas download
python -m pmp22_atlas analyze
python -m pmp22_atlas build --output build
python -m unittest discover -s tests -v
```

The build command validates every input before producing:

* `build/evidence_matrix.tsv`, one row per region/evidence layer;
* `build/candidate_panel.tsv`, the ranked, explicitly non-probabilistic panel;
* `build/source_overlap.dot`, a Graphviz source/exposure graph; and
* `build/atlas_summary.json`, machine-readable counts and unresolved gaps.

## Data contract

All genomic exports use zero-based, half-open intervals.  Unknown coordinates
are represented by an empty value, never invented.  `evidence_status` is one of
`measured`, `predicted`, `proposed`, or `missing`; those categories must not be
collapsed.  Candidate priority is an ordinal, rule-based triage score and must
not be interpreted as probability of efficacy.

The committed seed records capture only facts supported by the handoff.  Rows
with unresolved human source, donor, interval, or exposure metadata stay
explicitly unresolved in `data/open_gaps.md`.  Raw assay ingestion and native
human perturbation remain acquisition/experimental dependencies.

## Executed processed-data analysis

The pinned GEO peak/signal tracks, GEO SOFT records, FANTOM SDRF, and UCSC rn5
RefGene table have been downloaded and checksum-verified. Raw downloads are
excluded from Git because they total roughly 160 MB; their URLs and SHA-256
digests are committed in `data/download_manifest.tsv`. The reproducible locus
slice and results are committed under `data/derived/`.

The analysis uses the RefGene Pmp22 span plus 1 Mb on each side. It is a
descriptive overlap/signal calculation on submitter-processed tracks, not a
differential test. See `data/derived/analysis_report.md` for results and limits.

## Layout

* `data/` — versioned manifests and evidence tables.
* `candidate_dossiers/` — falsifiable hypotheses and decisive experiments.
* `pmp22_atlas/` — dependency-free validation and deterministic build CLI.
* `tests/` — schema, coordinate, foreign-key, and build regression tests.
* `benchmark/README.md` — frozen benchmark contract and stop conditions.
