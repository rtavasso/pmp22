# Project A — PMP22 regulatory atlas

The public-data computational release includes a human interval atlas,
donor-resolved CAGE, Schwann and bulk-nerve evidence, published perturbation
effects, an executed genomic benchmark, and seven experimental dossiers.
**Read [report.html](report.html)** for the findings and their limits.

The original project_a_regulatory_atlas.html remains the research contract.
The earlier branch contained a rat processed-track pilot. This continuation
does not assert native human regulatory causality, a human state effect or
therapeutic efficacy. [COMPLETION_PLAN.md](COMPLETION_PLAN.md) provides the
acceptance audit; [METHODS.md](METHODS.md) explains the analyses.

## Reproduce

Python 3.13.2 was used. Create a virtual environment and install
requirements-analysis.txt, then run these commands with its Python executable:

    python scripts/restore_resources.py
    python -m pmp22_atlas download
    python -m pmp22_atlas analyze
    python -m pmp22_atlas human
    python -m pmp22_atlas context
    python scripts/curate_evidence.py
    python scripts/build_dossiers.py
    python -m pmp22_atlas benchmark
    python -m pmp22_atlas supplemental
    python -m pmp22_atlas comparability
    python -m pmp22_atlas build --output data/release
    python -m pmp22_atlas report
    python -m unittest discover -s tests -v

For example, python -m venv .venv followed by the environment's Python
-m pip install -r requirements-analysis.txt installs the analysis packages.
On Windows that executable is .venv/Scripts/python.exe; on Unix it is
.venv/bin/python. uv pip install --python <executable> -r
requirements-analysis.txt is an alternative when pip is unavailable.

A fresh full replay downloads several GB of public reference/assay data.
Compressed public metadata snapshots preserve the acquired API responses
including changing timestamps. Large references, assays and papers stay outside
Git. resource_lock.json and download_manifest.tsv pin SHA-256; ENCODE downloads
also check publisher MD5. Changed remote bytes fail explicitly. The original
benchmark split and protocol refuse overwrite on a mismatching replay.

The report, measurements, predictions, split IDs, model coefficients,
effects and source records are committed and readable without raw downloads.
Report generation uses committed tables and does not refit models.

## Results

- Six historical human constructs map uniquely and reciprocally to GRCh38;
  nine human intervals include promoter families and the distal envelope.
- Distal C and the intronic element have stringent pooled Schwann ATAC support.
  A and B appear only at the broader Schwann peak threshold.
- Cultured donor1 has 22 P1 / 285 P2 CAGE tags; donor3 has 1 / 122. Donor2
  has only 0 / 1 tags and cannot support a precise promoter comparison.
- Eleven bulk-tibial assays represent four donors, overlapping the two ATAC
  donors. They are not eleven independent donor validations.
- The sequence composite reaches AP 0.966 on 2,120 held-out genomic windows
  and 0.784 on a different, 600-example GC-matched population. Matching retains
  28.3% of examples, changes promoter composition and leaves CpG differences;
  the score change cannot be attributed to GC alone. Distance-to-TSS reaches 0.771 there;
  the paired AP difference CI includes zero. Sequence-specific superiority,
  disease prediction and enhancer-to-gene links are not established.
- The matched subset beyond 2 kb from annotated TSSs has 139 positives and
  247 backgrounds: frozen sequence AP 0.580, AUROC 0.746, prevalence 0.360.
  This post-review diagnostic is not independent validation of enhancers.
- Jones 2011's intronic 50% effect is reduced EGR2 fold induction, each construct
  normalized to its own no-EGR2 baseline; it is not a directly comparable loss
  of induced output. See REVIEW_RESOLUTION.md for all adversarial corrections.
- Standard Borzoi splits place PMP22 and its centered 524,288 bp contexts
  in training fold7. No independent PMP22 test or foundation-model
  head-to-head comparison is claimed.

## Deliverables

data/human/ contains observations, coordinate audits, BED/FASTA, bulk context,
rodent orthology and QC. data/assay_effects.tsv holds 22 typed published
records including negative mutations. benchmark/results/ holds predictions,
uncertainty and controls. data/release/ contains the evidence matrix and
candidate panel. candidate_dossiers/ contains seven hypotheses and the B/D
handoff. data/open_gaps.md identifies the remaining scientific dependencies.

Exports are 0-based and half-open. Historical printed coordinates retain a
possible 1 bp convention ambiguity until inserts/junctions are confirmed.
Reference FASTA is not a sequenced plasmid. Candidate priority is ordinal,
never a probability or a predicted human effect size.
