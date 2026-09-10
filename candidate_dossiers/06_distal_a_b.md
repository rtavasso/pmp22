# Distal A versus B: discriminate factor-dependent regulation

- **Priority:** 3
- **Direction:** repression / mechanistic comparison
- **Decisive experiment:** Test A and B separately with matched local inhibition and factor-occupancy measurements.
- **Region:** PMP22_DISTAL_A; GRCh38 chr17:15409812-15410108 (0-based, half-open; gene strand minus)
- **Status:** Proposed experiment; no native human effect measured.

## Evidence and rationale

A and B are broad-peak-only Schwann intervals; both have bulk H3K27ac. A has EP300 overlap while B lacks the selected EP300 peak. Published S16 reporter effects differ: A SOX10-site mutation approximately 85% loss; B EGR2_1 approximately 45% loss. B SOX10/EGR2_2 mutations were retained negative results.

## Falsifiable hypothesis

A and B may have different native regulatory contributions depending on the confirmed SOX10/EGR2 context. Different reporter mutations and denominators do not establish their relative native strength.

## Decisive contrast and interpretation

Compare region-specific inhibition without first perturbing the whole transcription-factor program. Loss of reporter activity alone does not establish native human regulation; factor knockdown alone can change cell identity.

## Main uncertainty

B interval: GRCh38 chr17:15405970-15406367 (0-based, half-open; gene strand minus). Peak-threshold sensitivity is uncertainty, not inactivity. Published mutation sequences require supplement/plasmid confirmation before recreation.

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
