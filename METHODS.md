# Methods and interpretation

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

One training-label permutation is a sanity check, not a significance test.
Conservation, shifted windows and additional CIs were declared after the first
test and are labeled supplemental. Released labels were called genome-wide
by the source pipeline: this is a retrospective fixed-label benchmark, not
end-to-end processing with train-only peak/bias estimation. No count/profile,
foundation-model head-to-head, human-state or native-perturbation prediction
is claimed. Candidates are not ranked by the classifier.
