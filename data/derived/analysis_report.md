# Executed processed-data analysis

## Inputs and interval

All 12 inputs in `data/download_manifest.tsv` passed SHA-256 verification. UCSC
rn5 RefGene places Pmp22 at `chr10:49316999-49346994` (+ strand).
The frozen analysis window is `chr10:48316999-50346994`: the unioned
RefGene span with a 1 Mb flank on both sides.

## Results

Seven of eight processed tracks contain records in the window. The PNS SOX10
peak set contains **35**
overlapping peaks, versus **3**
in the CNS peak set. The differential peak files contribute
**21 PNS** and
**0 CNS**
records. The P15 sciatic-nerve H3K27ac set contributes
**28** peaks.

Across covered bases (zeros/unreported bases excluded), submitter bedGraph signal
averages **4.081578**
for PNS SOX10, **0.775705**
for CNS SOX10, and
**3.609413**
for PNS H3K27ac. These track scales are not assumed comparable.

The FANTOM SDRF extraction independently recovers three distinct Schwann records:
CNhs12073/donor1, CNhs12345/donor2, and CNhs12621/donor3. It does not resolve the
Borzoi label/path mismatch. Donor-specific CAGE measurements are now available
in the separate human analysis under `data/human/`.

## Interpretation boundary

The result reproduces descriptive processed-track support around rat Pmp22 and
confirms the FANTOM donor metadata. It does **not** estimate a state effect: the
bedGraphs are combined tracks, peak counts ignore width and threshold differences,
and the H3K27ac experiment has one IP plus input rather than two biological
replicates. It is not evidence of native human causality, therapeutic selectivity,
or a numerical PMP22 RNA/protein response.
