# Frozen analysis contract

## Primary contrast

For each existing native perturbation dataset that passes metadata/QC review,
estimate the assay's treatment-versus-matched-control effect on PMP22 RNA.  RNA,
promoter usage, protein, chromatin, neighbor genes, cellular identity, and
function are separate outcomes.  No suitable native human perturbation dataset
is currently present, so this contrast is **pending**, not zero.

## Exploratory contrasts

1. Rat sciatic nerve versus spinal cord SOX10 occupancy (descriptive if pools
   cannot establish biological replication).
2. Injury versus sham state evidence, conditional on recovering its design.
3. Human donor-resolved promoter usage and accessibility, conditional on source
   reconciliation and an estimable donor/study design.
4. Within-profile genomic prediction, clearly distinguished from donor/state
   generalization.

## Statistical commitments

* Freeze chromosome/region-block and donor/study splits before feature fitting.
* Define peaks and fit bias correction on training data only.
* Report effect sizes and uncertainty; never substitute inputs, strands, cells,
  or sequencing runs for biological replicates.
* Control multiplicity within each declared outcome family and report all tested
  candidates, including null and negative-control results.
* Report effective independent elements and studies beside nominal observations.
* Treat study-state confounding as non-identifiability, not a batch-correction
  problem. Restrict the conclusion when the design does not cross study/state.
* Benchmark distance, motif, conservation, and measured-accessibility baselines
  before deep models, using the same split and outcome.

## Precision and deviation policy

Numerical precision targets remain pending until pilot variance and independent
sample counts are known.  Power is not retroactively defined from observed
significance.  Every change to outcomes, exclusions, transforms, splits, or
models is dated in a deviation log before inspecting affected test results.

