# Distal C: test a within-donor state interaction

- **Priority:** 5
- **Direction:** state-conditional regulation
- **Decisive experiment:** Cross mature/myelinating and repair-like states with C perturbation within each donor.
- **Region:** PMP22_DISTAL_C; GRCh38 chr17:15377645-15378054 (0-based, half-open; gene strand minus)
- **Status:** Proposed experiment; no native human effect measured.

## Evidence and rationale

C is a human accessibility-supported interval chosen for a falsifiable state-interaction test. Current C-specific differential evidence is unavailable: two partially mapped rat sham peaks overlap 63 and 127 bp of its human interval, but no submitted differential peak overlaps exact C. Whole-nerve peak counts across a 2 Mb neighborhood do not increase its priority or establish regulation within Schwann cells.

## Falsifiable hypothesis

The effect of C inhibition on PMP22 differs between experimentally confirmed human states.

## Decisive contrast and interpretation

Estimate the intervention-by-state interaction, with both states represented in every donor and balanced batches. Compare targeting occupancy across states. A study-versus-state comparison cannot identify this interaction.

## Main uncertainty

Human CAGE and adult ATAC cannot establish a state trajectory. Whole-nerve injury chromatin also confounds regulation within cells with cell composition. Orthology counts are conditional on mapping thresholds; inspect filtered candidates, not only accepted counts.

## Shared measurement contract

Use independently genotyped human Schwann donor lines with the state confirmed
in each differentiation. Cross donors with intervention and state; do not assign
one study or donor exclusively to a state. Include non-targeting controls, an
effector-only control, the proposed matched DNA intervals, and at least two
independent targeting reagents per region. Reagents and on-target occupancy
need validation before interpreting a negative expression result.
Include the two proposed accessible comparison loci alongside the closed DNA
controls. Confirm comparable recruitment and local chromatin modulation in the
actual donor/state; absence of recruitment at a closed site is not a functional
negative control for inhibition at an accessible enhancer. Accessible comparison
loci are also untested and may regulate other genes or PMP22.

Primary endpoint: donor-paired change in total PMP22 RNA relative to the matched
control. Measure P1/exon1A and P2/exon1B initiation separately; measure protein
as its own later endpoint. Keep early local neighbor effects separate from later
downstream rescue responses. Monitor human TVP23C, CDRT4, TEKT3 and other genes
in the versioned locus table, plus SOX10/EGR2, MPZ/MAG and viability. The rodent
Tvp23b control must not silently be relabeled as human TVP23C.

Use donor as the biological unit and differentiation batch as a nested factor.
Technical wells and guides are not donor replicates. Start with a feasibility
pilot spanning at least three donors and two differentiation batches; this is
not a powered efficacy study. Freeze the primary comparison and minimum useful
RNA effect before collecting the confirmatory data. Use pilot between-donor
variance to choose the confirmatory donor count for a target CI width. Report
all candidates, null outcomes and exclusions; correct the confirmatory family
of candidate RNA tests for multiplicity. No magnitude or probability of human
response is assigned by this atlas or its accessibility classifier.

A local PMP22 effect preceding loss of differentiation supports the regulatory
hypothesis. A generalized identity/toxicity effect, failure of targeting, or
an imprecise null requires a different interpretation. A selective early
effect can still produce a broad downstream disease-response program.


## Handoff

Project D receives this interval, the BED/reference FASTA, source assay effects,
negative controls, donor/state design and endpoint contract. Guide sequences,
on-target binding, genotyping and delivery need experimental design and review.
Project B receives promoter/UTR and transcript-annotation ambiguity. Printed
historical coordinates retain a possible 1 bp convention uncertainty; reference
FASTA is not a sequenced plasmid insert.
