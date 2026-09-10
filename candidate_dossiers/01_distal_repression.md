# Distal C: test native human repression

- **Priority:** 1
- **Direction:** repression
- **Decisive experiment:** Perturb distal C and measure early total/P1/P2 RNA against matched controls.
- **Region:** PMP22_DISTAL_C; GRCh38 chr17:15377645-15378054 (0-based, half-open; gene strand minus)
- **Status:** Proposed experiment; no native human effect measured.

## Evidence and rationale

C has one stringent and one broad Schwann ATAC interval, bulk H3K27ac and EP300 support, and a published human-sequence reporter with an approximately 65% loss after SOX10-site mutation in rat S16 cells. Those observations support choosing the interval; 65% is not a human repression forecast.

## Falsifiable hypothesis

Local inhibition of C reduces native PMP22 before changing Schwann identity.

## Decisive contrast and interpretation

An unchanged PMP22 endpoint with verified local inhibition and adequate precision would argue that C is dispensable in the tested human state. Joint C plus intronic inhibition can distinguish redundancy from an inactive C.

## Main uncertainty

Orthologous rat sham peaks overlap C, but no submitted sham-enriched differential peak overlaps its exact reporter interval. A mouse deletion covers a much larger, partly conserved domain.

## Shared measurement contract

Use independently genotyped human Schwann donor lines with the state confirmed
in each differentiation. Cross donors with intervention and state; do not assign
one study or donor exclusively to a state. Include non-targeting controls, an
effector-only control, the proposed matched DNA intervals, and at least two
independent targeting reagents per region. Reagents and on-target occupancy
need validation before interpreting a negative expression result.

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
