# Executed benchmark

protocol.json freezes version 1, seed 20260910. splits.tsv records each
example, label, chromosome, sequence hash, GC, TSS distance and split.
results/dataset_lock.json pins the split and protocol. Parameters are JSON,
not executable pickle files.

The measurement is binary accessibility in one released pooled Schwann profile:
an exploratory retrospective genomic task, with chromosome holdouts but no
donor holdout. Source-pipeline labels are fixed external calls, not train-only
raw-read peak estimates. See amendments.md for the declared deviation and
post-test supplemental diagnostics.

Input tiers: sequence GC/CpG, sequence motif maxima, sequence composite
(GC/CpG + motifs + canonical 4-mers), annotation-based TSS distance, and
supplemental multispecies conservation. Accessibility is the outcome, so
using measured accessibility to predict itself would leak the answer.
No matched RNA/perturbation test supports an accessibility-input comparator.

Scaling and logistic fitting use chr1/2; C uses chr3; chr17/22 are test-only.
AP is primary; AUROC, Brier and calibration are secondary. The 50:50 sample
is not natural genomic prevalence. All models use the same held-out IDs.

The sequence composite reaches AP 0.966 (CI 0.955-0.976), versus GC/CpG
0.944. In the GC-matched set, sequence reaches 0.784 (0.732-0.833), distance
0.771 (0.727-0.813), GC/CpG 0.648, motifs 0.674 and conservation 0.689.
The sequence-minus-distance CI spans zero. The permutation control reaches
0.473 overall and 0.509 when matched. 500 paired 1 Mb block draws estimate
uncertainty conditional on the single source.

Standard Borzoi replicas share test fold3/validation fold4. PMP22 and centered
524,288 bp inputs overlap training fold7; these checkpoints cannot certify
an independent PMP22 test. No model is selected by PMP22 scores, and candidates
are not classifier-ranked. AlphaGenome, ChromBPNet, CAGE-profile and mutation
benchmarks were not executed: their matched-data/independence gates remain
unmet rather than being replaced with fabricated or incompatible scores.
