"""Write the computational release status and reproducibility contract."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import ROOT
from pmp22_atlas.human import table

def main():
    docs={
'README.md':"""# Project A — PMP22 regulatory atlas

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
""",
'METHODS.md':"""# Methods and interpretation

## Human reference and mapping

UCSC hg38 primary chromosomes and GENCODE50 tables are byte-pinned. No patch
contigs were used. Jones 2011/2012 printed hg18 reporter coordinates are
interpreted as 1-based closed, subtracting one from the start for BED.
All six intervals have full chain coverage, a unique target and exact
reciprocal endpoint mapping. Actual cloned inserts remain unverified, with
a possible 1 bp convention uncertainty rather than asserted plasmid identity.

PMP22 is on the minus strand. P1/exon1A and P2/exon1B use the last genomic
exon of coding GENCODE50 transcript families. Gene-oriented reference FASTA is
reverse-complemented. The distal envelope spans three human reporters; it is
not a tested human deletion. Mouse SE mapping covers 54.5% of its source span.
Its gap-containing target envelope is not interchangeable with aligned blocks.
Rat injury overlap uses aligned blocks, not the mapped envelope.
rat_mapping_overlap_audit.tsv retains low-coverage and ambiguous candidates;
region_rat_injury_support.tsv exposes filtered counts and evaluation status.
A zero accepted count is not a complete census of biological activity.

## CAGE and annotation sensitivity

Each FANTOM hg19 CTSS file is counted separately. All tag weights form the
library denominator. Count minus-strand tags in each first-exon-family window,
uniquely mapped from hg38 to hg19. The primary window adds 100 bp on each side;
upstream-flank sensitivities add 250 or 500 bp. TPM means CTSS tags per million
total tag weights, not transcript TPM or mature-RNA isoform abundance.

All three GENCODE window sizes agree. RefSeq coding first-exon boundaries
give P2 donor1=284 rather than 285 tags; all other counts agree. Cultured
samples do not certify myelinating identity. Donor2 has 232,185 library tags
and only one promoter tag: do not infer that P1 is absent.

## Human chromatin, source QC and bulk RNA

ENCFF602YVR supplies IDR-thresholded peaks; ENCFF632FMT supplies broader
pseudoreplicated peaks. Deduplicate exact chromosome/start/end coordinates.
These thresholds are sensitivity views of one pooled product, not replicates.
The IDR file is a pseudoreplicate comparison, not between-donor replication.
The flattened matrix carries source_interval_ids, evidence_family_id and
donor_ids. For example, one 25,903 bp H3K27ac interval supplies five region
overlaps. Families are dependency labels, not independent replicate counts.
Source alignments trace to ENCDO793LXB (female, 53) and ENCDO271OUW (female, 51).
The portal adult/child label does not identify an additional child donor.
Whole-library mapping/complexity metrics cannot substitute for Schwann-subset
FRiP, cell count or TSS enrichment; the peak-file QC lists are empty.

Portal audit flags are retained. Matching-MD5 flags resolve to an alternate
leg-Schwann annotation ENCSR221DFU, which contains identical peak bytes.
It is duplicate data, not replication. CATlas kai199 equivalence and the Borzoi
FANTOM label/path mismatch remain unresolved and cannot certify independence.

For all eleven bulk-tibial leads, choose a released GRCh38 BED/TSV by output
class and newest annotation/release before examining the locus. Selection,
donor, library, biosample and control records are exported. Four RNA assays
use GENCODE29 gene quantifications. Retain source gene IDs and match symbols
through current transcript IDs. Total PMP22 uses ENSG00000109099. Historical
ENSG00000230971 shares a transcript now assigned to PMP22; retain it separately
instead of silently combining models. RSEM posterior intervals do not measure
between-donor biological variance. RNA TPM is not compared to CAGE TPM.
Bulk chromatin is not purified Schwann evidence; CTCF is not a measured loop.

## Rodent data and published effects

The original twelve pilot inputs are checksum reverified and reanalyzed in
the rn5 gene span plus/minus 1 Mb. Track counts and covered-base signal are
descriptive; different scales and peak widths are not assumed comparable.
GSE63103 has two H3K27ac IP pools per sham/injury condition plus inputs.
Reuse author differential labels instead of refitting a differential test to
pooled intervals. Eight precision-ambiguous scientific-notation coordinates
are rejected, all outside the locus window. Submitted 99-ended BED coordinates
retain an unresolved source convention.
Whole-nerve injury measurements leave cell composition unresolved. The 83/38
window peak counts and partial sham-only C overlaps do not identify a C-specific
state effect. No submitted differential peak overlaps the exact C reporter.

ENA runs are summed within GEO libraries, not counted as biological replicates.
PNS SOX10 libraries have 29.6 versus 154.6 million raw reads. Archive read depth
is not mapped usable depth. Animals per pool remain unresolved in deposited
metadata. GSE64971 is one IP plus input, not n=2.

Effects are author-text approximations, exact counts or qualitative directions.
No figure is digitized. Numerical SDs, native mouse magnitudes and individual
values remain blank when unavailable. Negative mutations and unchanged controls
are retained; nonsignificance is not equivalence. The rat S16 deletion line has
approximately three chr10 copies. Three deletion clones are compared to five
controls; its unaffected reporter allele is an internal control. Human
whole-PMP22-copy deletion is excluded from noncoding-only effect labels.
Jones 2011 Figure 3 measures EGR2 fold induction: each construct's +EGR2
activity divided by its own no-EGR2 baseline. AE011 is an approximate 50%
reduction in that ratio, not induced output. Figure 3 reports n=6 but does not
resolve measurement independence. Mouse RNA fields retain author significance
without asserting an ANOVA-only test absent from the panel's documentation.
Different reporter mutations, induction ratios and native RNA effects are not
pooled into a potency ranking. Source pairing and conduction-endpoint wording
remain unresolved; no new P values or loss percentages are reconstructed.

## Benchmark

The frozen protocol, split and dated amendments are in benchmark/. The outcome
is peak-containing versus sampled-background 600 bp windows in one fixed pooled
profile, not validated enhancer function. Chr1/2 train, chr3 validation and
chr17/22 test contain 5,916/2,166/2,120 examples, balanced within chromosome.
Test examples occupy 119 1 Mb blocks in one source, not 2,120 donors.

Overlap clustering, blacklist/N filtering, a broad-peak exclusion buffer and
reverse-complement sequence-hash deduplication precede fitting. Standardization
and fitting use training examples only; regularization uses validation AP.
The sequence composite contains GC, CpG, three motif maxima and canonical
4-mer frequencies. It has no TSS distance and is not a pure 4-mer ablation.
Distance uses gene annotation; conservation uses a multispecies alignment.
Those input tiers are distinct.

AP is primary; AUROC, Brier score and calibration are secondary. AP=0.5 is the
constructed 50:50 prevalence, not genome prevalence. Scores are not calibrated
probabilities across arbitrary loci or probabilities of experimental success.
500 paired genomic-block bootstrap draws condition on this source and sampling.
The GC-matched subset has 600 examples in 0.02 GC bins within chromosome,
occupying 110 blocks. It is evaluated without refitting. Sequence-minus-distance
AP is 0.013 with CI -0.036 to 0.062, so superiority is not established.
Matching retains 28.3% of examples and changes the fraction of positives within
2 kb of an annotated TSS from 78.3% to 53.7%. It does not match CpG frequency.
The comparability command recomputes CpG from checksum-verified reference
sequences and scores frozen predictions in explicit strata. Original versus
matched AP is not a causal decomposition of GC effects. Subgroup AP must carry
its own prevalence; the matched beyond-2kb stratum has prevalence 0.360 and
sequence AP 0.580/AUROC 0.746. These post-review analyses do not change v1.

One training-label permutation is a sanity check, not a significance test.
Conservation, shifted windows and additional CIs were declared after the first
test and are labeled supplemental. Released labels were called genome-wide
by the source pipeline: this is a retrospective fixed-label benchmark, not
end-to-end processing with train-only peak/bias estimation. No count/profile,
foundation-model head-to-head, human-state or native-perturbation prediction
is claimed. Candidates are not ranked by the classifier.

## Experimental controls

Two closed, GC/length-matched intervals are retained for their distinct role.
Two accessible comparison loci now use 409 bp summit-centered strict ATAC
windows, GC difference <=0.03, outside the PMP22 +/-1 Mb neighborhood,
GENCODE50 basic transcripts and 2 kb TSS windows, with blacklist/N filtering.
Selection is deterministic by GC difference then distance, with >=5 kb spacing.
All four DNA sites remain untested. Absence of a known PMP22 link is not proof
of no effect. Confirm recruitment and local modulation at every comparison
site in the experimental donor/state; closed controls alone do not match an
accessible enhancer's recruitment context. Add non-targeting/effector controls.
""",
'benchmark/README.md':"""# Executed benchmark

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
""",
'data/open_gaps.md':"""# Remaining scientific evidence

1. Native human noncoding perturbation on RNA, protein, neighbors and function.
   Seven dossiers specify the new experiments. This is the main scientific
   gate, not a missing software command.
2. Replicated human states crossed within donors. Cultured CAGE versus adult
   tissue ATAC cannot isolate state; fetal CATlas labels do not solve the design.
3. Matched human Schwann contact/target evidence. Bulk CTCF, EP300 and distance
   remain priors, not measured loops.
4. Historical cloned inserts and rat deletion junctions. Printed human
   intervals map fully but retain possible 1 bp convention ambiguity. Mouse
   SE mapping covers only 54.5% of its span.
5. Original numerical assay uncertainties and mouse effect magnitudes.
   Author approximations are retained; no figure values or meta-analysis
   precision are invented.
6. CATlas kai199/fetal-source reconciliation and Borzoi FANTOM donor-label
   recipe. PMP22 training exposure itself is resolved. Independent profiles
   or perturbation labels are needed for broader model claims.
7. Donor-specific Schwann subset cell counts, FRiP, TSS enrichment and bias
   estimates. Original source QC is preserved but raw-read processing and
   count/profile modeling were not performed.
8. Novelty and confirmatory precision. Known enhancers are not rediscovered
   as new regulators. Set a minimum useful effect and use pilot donor variance
   before choosing a powered confirmatory sample size.

Resolved since the pilot: human peak access and ATAC donor paths; three CAGE
measurements; historical human mappings; GSE63103 injury files; rat checksum
replay; eleven bulk leads/shared donors; a benchmark with controls and uncertainty.
""",
'data/permissions.md':"""# Public access and attribution

Resources were served publicly without controlled-access credentials:
GEO/ENA processed files and metadata, ENCODE releases, FANTOM, UCSC,
JASPAR and official model source inventories. Original URLs and timestamps
are retained in the resource lock and download manifest. Observe provider
attribution and applicable reuse terms:

- https://www.encodeproject.org/help/data-use-policy/
- https://fantom.gsc.riken.jp/5/
- https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html
- https://genome.ucsc.edu/license/
- https://jaspar.elixir.no/

No controlled-access GTEx genotypes were acquired. Published PDFs are not
redistributed in Git; extracted facts have source locators. Public metadata
snapshots preserve provenance, not a new license for source content.
""",
'analysis_contract.md':"""# Analysis contract

The primary scientific contrast remains intervention versus matched control
for native PMP22 RNA, separately from promoter usage, protein, chromatin,
neighbors, identity and function. Published rodent perturbations support
their own context. No new native human noncoding effect is estimated.

Executed observational analyses include cultured CAGE, threshold-sensitive
Schwann ATAC, bulk-nerve context, and rat tissue/injury maps. METHODS.md and
sample/source manifests define their units and limits. The benchmark protocol
and dated amendments distinguish its genomic classification task from raw-data
profile modeling, state effects and experimental response.

Source/state confounding is non-identifiability, not solved by batch adjustment.
Genomic blocks are not biological donors. Dossiers specify donor-paired,
crossed-state experiments, local/identity controls and a precision-design gate.
An imprecise null is not equivalence and retrospective significance is not
power. The original HTML remains the broader acceptance contract.
""",
    }
    for name,text in docs.items():(ROOT/name).write_text(text,encoding='utf-8')
    gates=[
        ('human_intervals','completed_with_caveat','Six mapped constructs; nine human intervals; BED/FASTA and reciprocal audit','Cloned inserts/junctions need confirmation'),
        ('human_measurements','completed_descriptively','Three CAGE donors; Schwann ATAC; eleven bulk assays/four donors','Pooled/cross-study data do not identify human states'),
        ('source_independence','completed_for_analyzed_sources','ATAC donors; bulk overlap; duplicate annotation; GEO runs','CATlas aliases and FANTOM model recipe unresolved'),
        ('published_effects','completed_with_missing_values_flagged','22 typed records including approximate values and negative mutations','Original numerical uncertainties and human noncoding effects absent'),
        ('rodent_state_context','completed_descriptively','Twelve inputs replayed; injury files/QC/orthology','No new differential or transported human effect'),
        ('genomic_benchmark','completed_for_binary_task','Frozen split; baselines; conservation; uncertainty and controls','No count/profile or independent human-effect benchmark'),
        ('Borzoi_eligibility','resolved_for_standard_checkpoints','PMP22 and centered full input contexts are in training fold7','No independent PMP22 test from these models'),
        ('candidate_handoff','completed','Seven dossiers; two closed and two accessible untested DNA controls; B/D handoff','Lab targeting, genotyping, recruitment/local modulation, power and outcomes remain'),
        ('native_human_causality','requires_new_evidence','Not asserted','Matched noncoding perturbation assays'),
        ('human_state_generalization','requires_new_evidence','Not asserted','Crossed donor/state design'),
        ('human_contact_assignment','unresolved','No measured loop asserted','Matched contact or perturbation evidence'),
        ('publication_novelty_and_precision','not_certified','Integration and benchmark limitations documented','Independent contribution and precision assessment'),
    ]
    rows=[dict(gate=a,status=b,delivered=c,remaining=d) for a,b,c,d in gates]
    table(ROOT/'data/acceptance_audit.tsv',rows)
    text="# Project A acceptance audit — 10 September 2026\n\n"
    text+="The computational public-data release is delivered. The original scientific contract is not fully satisfied: new native human evidence and replicated human-state validation remain required. Software checks cannot close those gates.\n\n"
    text+="| Gate | Status | Delivered | Remaining |\n|---|---|---|---|\n"
    for r in rows:text+='| '+' | '.join(r.values())+' |\n'
    (ROOT/'COMPLETION_PLAN.md').write_text(text,encoding='utf-8')
    checks=dict(checks=[
        dict(id='schema',command='python -m pmp22_atlas validate',scope='manifest integrity'),
        dict(id='tests',command='python -m unittest discover -s tests -v',scope='coordinates, sources, splits and results'),
        dict(id='processed_replay',status='executed',scope='public processed inputs; not read-level validation'),
        dict(id='genomic_benchmark',status='executed',scope='one-profile binary classification'),
        dict(id='native_human_perturbation',status='not_performed',scope='new causal evidence required')],
        boundary='Scientific acceptance is separate from software checks; see data/acceptance_audit.tsv')
    (ROOT/'audit_checks.json').write_text(json.dumps(checks,indent=2)+'\n')

if __name__=='__main__':main()
