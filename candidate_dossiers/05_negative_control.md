# Closed and accessible DNA controls: test distinct assumptions

- **Priority:** 7
- **Direction:** control
- **Decisive experiment:** Compare closed and accessible DNA sites alongside non-targeting and known regulatory controls, verifying recruitment and local modulation.
- **Region:** PMP22_DISTAL_C; GRCh38 chr17:15377645-15378054 (0-based, half-open; gene strand minus)
- **Status:** Proposed experiment; no native human effect measured.

## Evidence and rationale

Two 409 bp intervals were selected without expression outcomes, with GC within 0.03 of C, outside broad ATAC plus a 1 kb buffer, at least 10 kb from known reporter/promoter intervals, outside annotated transcripts and the ENCODE blacklist: NEG_C_1 chr17:15354300-15354709; NEG_C_2 chr17:15426800-15427209.

## Falsifiable hypothesis

Two additional accessible comparison sites overlap stringent pooled Schwann ATAC, match C length and GC within 0.03, and lie outside the PMP22 neighborhood, annotated transcripts and 2 kb TSS windows: ACCESSIBLE_C_1 chr17:62694510-62694919; ACCESSIBLE_C_2 chr17:34199943-34200352. Targeting each candidate control does not change native PMP22 under the same assay conditions; this is a hypothesis, not an established negative.

## Decisive contrast and interpretation

Define an equivalence margin and collect enough donor-level precision before calling a control functionally negative. A detected effect rejects the negative assumption and requires local-gene/identity investigation.

## Main uncertainty

Closed sites can fail recruitment; accessible sites may have regulatory targets. Confirm recruitment and effector activity for every site. These controls test different assumptions and are not benchmark-validated negatives.

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
