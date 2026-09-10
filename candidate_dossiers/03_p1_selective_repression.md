# P1: test promoter selectivity in the right state

- **Priority:** 4
- **Direction:** promoter-specific repression
- **Decisive experiment:** Compare P1-directed inhibition with P2-directed and non-targeting controls in the same donor/state.
- **Region:** PMP22_P1; GRCh38 chr17:15265153-15265360 (0-based, half-open; gene strand minus)
- **Status:** Proposed experiment; no native human effect measured.

## Evidence and rationale

P1 is an annotated promoter-associated first exon but has no overlapping released Schwann peak at either threshold. Cultured CAGE has only 22/1 P1 tags in the two informative donors versus 285/122 P2 tags.

## Falsifiable hypothesis

P1 contributes more strongly in a confirmed myelinating state and can be modulated without the same P2 effect.

## Decisive contrast and interpretation

First establish measurable P1 initiation. A null in P2-dominant culture with little P1 activity does not test the myelinating-state hypothesis. Direct promoter targeting is a mechanistic control, not proof of enhancer selectivity.

## Main uncertainty

GENCODE50 and RefSeq promoter boundaries differ; counts are robust here but the exact targeting window requires TSS confirmation.

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
