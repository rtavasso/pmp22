"""Compile seven evidence-linked experimental hypotheses from the measured atlas."""
import csv
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import ROOT

def read(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def main():
    regions={r['region_id']:r for r in read(ROOT/'data/regions.tsv')}
    def interval(rid):
        r=regions[rid];return f"GRCh38 {r['chrom']}:{r['start']}-{r['end']} (0-based, half-open; gene strand minus)"
    negatives=read(ROOT/'data/human/candidate_negative_controls.tsv')
    accessible=read(ROOT/'data/human/candidate_accessible_controls.tsv')
    common="""
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
"""
    records=[
        ('01_distal_repression','Distal C: test native human repression',1,'repression','PMP22_DISTAL_C',
         'Perturb distal C and measure early total/P1/P2 RNA against matched controls.',
         'C has one stringent and one broad Schwann ATAC interval, bulk H3K27ac and EP300 support, and a published human-sequence reporter with an approximately 65% loss after SOX10-site mutation in rat S16 cells. Those observations support choosing the interval; 65% is not a human repression forecast.',
         'Local inhibition of C reduces native PMP22 before changing Schwann identity.',
         'An unchanged PMP22 endpoint with verified local inhibition and adequate precision would argue that C is dispensable in the tested human state. Joint C plus intronic inhibition can distinguish redundancy from an inactive C.',
         'Orthologous rat sham peaks overlap C, but no submitted sham-enriched differential peak overlaps its exact reporter interval. A mouse deletion covers a much larger, partly conserved domain.'),
        ('02_intronic_repression','Intronic enhancer: test a second regulatory route',2,'repression','PMP22_INTRONIC',
         'Test the intronic interval alone and jointly with distal C; distinguish transcription from splicing.',
         'The 243 bp reference interval maps reciprocally from the Jones 2011 reporter and overlaps stringent Schwann accessibility plus bulk H3K27ac/EP300. EGR2 site-4 mutation reduced EGR2 fold induction by approximately 50% in mouse B16/F10 cells. Each construct was normalized to its own no-EGR2 baseline; this is not a 50% reduction in induced output. Figure 3 reports n=6, with measurement independence unresolved.',
         'The intronic element contributes to native PMP22 transcription independently of distal C.',
         'Measure nascent initiation, common-exon RNA, splice junctions and protein. RNA loss driven by abnormal splicing or broad identity change would not establish enhancer-selective repression.',
         'Intronic targeting may alter elongation/splicing or spread repression. The larger transgenic construct is a separate interval and must not replace the small reporter without an explicit comparison.'),
        ('03_p1_selective_repression','P1: test promoter selectivity in the right state',4,'promoter-specific repression','PMP22_P1',
         'Compare P1-directed inhibition with P2-directed and non-targeting controls in the same donor/state.',
         'P1 is an annotated promoter-associated first exon but has no overlapping released Schwann peak at either threshold. Cultured CAGE has only 22/1 P1 tags in the two informative donors versus 285/122 P2 tags.',
         'P1 contributes more strongly in a confirmed myelinating state and can be modulated without the same P2 effect.',
         'First establish measurable P1 initiation. A null in P2-dominant culture with little P1 activity does not test the myelinating-state hypothesis. Direct promoter targeting is a mechanistic control, not proof of enhancer selectivity.',
         'GENCODE50 and RefSeq promoter boundaries differ; counts are robust here but the exact targeting window requires TSS confirmation.'),
        ('04_state_specific_distal','Distal C: test a within-donor state interaction',5,'state-conditional regulation','PMP22_DISTAL_C',
         'Cross mature/myelinating and repair-like states with C perturbation within each donor.',
         'C is a human accessibility-supported interval chosen for a falsifiable state-interaction test. Current C-specific differential evidence is unavailable: two partially mapped rat sham peaks overlap 63 and 127 bp of its human interval, but no submitted differential peak overlaps exact C. Whole-nerve peak counts across a 2 Mb neighborhood do not increase its priority or establish regulation within Schwann cells.',
         'The effect of C inhibition on PMP22 differs between experimentally confirmed human states.',
         'Estimate the intervention-by-state interaction, with both states represented in every donor and balanced batches. Compare targeting occupancy across states. A study-versus-state comparison cannot identify this interaction.',
         'Human CAGE and adult ATAC cannot establish a state trajectory. Whole-nerve injury chromatin also confounds regulation within cells with cell composition. Orthology counts are conditional on mapping thresholds; inspect filtered candidates, not only accepted counts.'),
        ('05_negative_control','Closed and accessible DNA controls: test distinct assumptions',7,'control','PMP22_DISTAL_C',
         'Compare closed and accessible DNA sites alongside non-targeting and known regulatory controls, verifying recruitment and local modulation.',
         'Two 409 bp intervals were selected without expression outcomes, with GC within 0.03 of C, outside broad ATAC plus a 1 kb buffer, at least 10 kb from known reporter/promoter intervals, outside annotated transcripts and the ENCODE blacklist: '+ '; '.join(f"{r['control_id']} chr17:{r['start']}-{r['end']}" for r in negatives)+'.',
         'Two additional accessible comparison sites overlap stringent pooled Schwann ATAC, match C length and GC within 0.03, and lie outside the PMP22 neighborhood, annotated transcripts and 2 kb TSS windows: '+ '; '.join(f"{r['control_id']} chr17:{r['start']}-{r['end']}" for r in accessible)+'. Targeting each candidate control does not change native PMP22 under the same assay conditions; this is a hypothesis, not an established negative.',
         'Define an equivalence margin and collect enough donor-level precision before calling a control functionally negative. A detected effect rejects the negative assumption and requires local-gene/identity investigation.',
         'Closed sites can fail recruitment; accessible sites may have regulatory targets. Confirm recruitment and effector activity for every site. These controls test different assumptions and are not benchmark-validated negatives.'),
        ('06_distal_a_b','Distal A versus B: discriminate factor-dependent regulation',3,'repression / mechanistic comparison','PMP22_DISTAL_A',
         'Test A and B separately with matched local inhibition and factor-occupancy measurements.',
         'A and B are broad-peak-only Schwann intervals; both have bulk H3K27ac. A has EP300 overlap while B lacks the selected EP300 peak. Published S16 reporter effects differ: A SOX10-site mutation approximately 85% loss; B EGR2_1 approximately 45% loss. B SOX10/EGR2_2 mutations were retained negative results.',
         'A and B may have different native regulatory contributions depending on the confirmed SOX10/EGR2 context. Different reporter mutations and denominators do not establish their relative native strength.',
         'Compare region-specific inhibition without first perturbing the whole transcription-factor program. Loss of reporter activity alone does not establish native human regulation; factor knockdown alone can change cell identity.',
         'B interval: '+interval('PMP22_DISTAL_B')+'. Peak-threshold sensitivity is uncertainty, not inactivity. Published mutation sequences require supplement/plasmid confirmation before recreation.'),
        ('07_deletion_hnpp_activation','Deletion-HNPP: test activation of the retained normal locus',6,'activation; deletion-HNPP only','PMP22_INTRONIC',
         'In deletion-HNPP lines with an intact retained coding allele, test graded local activation and measure RNA, protein and function.',
         'Human intronic and distal C accessibility provides a tractable regulatory hypothesis. Rodent SE deletions establish that reducing this regulatory system can lower Pmp22 and increase vulnerability, but reversing a deletion experiment does not establish activation rescue.',
         'Activating the retained normal locus can increase PMP22 while preserving differentiation and neighbor expression.',
         'Verify copy number and retained coding sequence, compare isogenic or donor-matched controls, and measure early transcription plus later protein and functional outcomes. Stop escalation if expression overshoots the predeclared healthy-reference interval or identity/viability deteriorates.',
         'Do not pool coding-variant HNPP with deletion-HNPP: nonselective activation could also raise a pathogenic allele. No patient-specific dose or therapeutic efficacy is inferred.'),
    ]
    for filename,title,priority,direction,rid,experiment,evidence,hypothesis,falsification,uncertainty in records:
        text=f"""# {title}

- **Priority:** {priority}
- **Direction:** {direction}
- **Decisive experiment:** {experiment}
- **Region:** {rid}; {interval(rid)}
- **Status:** Proposed experiment; no native human effect measured.

## Evidence and rationale

{evidence}

## Falsifiable hypothesis

{hypothesis}

## Decisive contrast and interpretation

{falsification}

## Main uncertainty

{uncertainty}

## Shared measurement contract
{common}

## Handoff

Project D receives this interval, the BED/reference FASTA, source assay effects,
negative controls, donor/state design and endpoint contract. Guide sequences,
on-target binding, genotyping and delivery need experimental design and review.
Project B receives promoter/UTR and transcript-annotation ambiguity. Printed
historical coordinates retain a possible 1 bp convention uncertainty; reference
FASTA is not a sequenced plasmid insert.
"""
        (ROOT/'candidate_dossiers'/(filename+'.md')).write_text(text,encoding='utf-8')
    print('Compiled seven candidate dossiers.')

if __name__=='__main__':main()
