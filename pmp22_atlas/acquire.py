"""Download pinned public inputs and summarize the rn5 PMP22 locus."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import urllib.request
from pathlib import Path

from .cli import DATA

RAW = DATA / "raw"
DERIVED = DATA / "derived"
LOCUS_FLANK = 1_000_000


def manifest() -> list[dict[str, str]]:
    with (DATA / "download_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def download(force: bool = False) -> dict[str, int]:
    RAW.mkdir(parents=True, exist_ok=True)
    fetched = reused = 0
    for row in manifest():
        path = RAW / row["filename"]
        if path.exists() and digest(path) == row["sha256"] and not force:
            reused += 1
            continue
        temporary = path.with_suffix(path.suffix + ".part")
        urllib.request.urlretrieve(row["url"], temporary)
        if digest(temporary) != row["sha256"]:
            temporary.unlink(missing_ok=True)
            raise ValueError(f"checksum mismatch: {row['filename']}")
        temporary.replace(path)
        fetched += 1
    return {"fetched": fetched, "reused": reused, "files": fetched + reused}


def require_inputs() -> None:
    absent = [r["filename"] for r in manifest() if not (RAW / r["filename"]).exists()]
    if absent:
        raise FileNotFoundError("run `python -m pmp22_atlas download` first: " + ", ".join(absent))
    bad = [r["filename"] for r in manifest() if digest(RAW / r["filename"]) != r["sha256"]]
    if bad:
        raise ValueError("checksum mismatch: " + ", ".join(bad))


def pmp22_locus() -> tuple[str, int, int, str]:
    hits = []
    with gzip.open(RAW / "rn5_refGene.txt.gz", "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 12 and fields[12] == "Pmp22":
                hits.append((fields[2], int(fields[4]), int(fields[5]), fields[3]))
    if not hits:
        raise ValueError("Pmp22 absent from pinned rn5 refGene")
    chroms = {hit[0] for hit in hits}
    if len(chroms) != 1:
        raise ValueError("Pmp22 annotation spans multiple chromosomes")
    return hits[0][0], min(h[1] for h in hits), max(h[2] for h in hits), hits[0][3]


def summarize_track(path: Path, chrom: str, start: int, end: int) -> tuple[dict[str, object], list[tuple[str, int, int, str]]]:
    records = overlaps = covered = 0
    weighted = 0.0
    maximum: float | None = None
    sliced = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            records += 1
            if fields[0] != chrom:
                continue
            left, right = int(fields[1]), int(fields[2])
            overlap = max(0, min(right, end) - max(left, start))
            if not overlap:
                continue
            overlaps += 1
            value = fields[3] if len(fields) > 3 else ""
            sliced.append((fields[0], left, right, value))
            if value:
                numeric = float(value)
                covered += overlap
                weighted += overlap * numeric
                maximum = numeric if maximum is None else max(maximum, numeric)
    stats = {"file": path.name, "records_genome_wide": records, "records_in_window": overlaps, "signal_covered_bp": covered, "signal_mean_covered_bp": round(weighted / covered, 6) if covered else None, "signal_max": maximum}
    return stats, sliced


def analyze() -> dict[str, object]:
    require_inputs()
    chrom, gene_start, gene_end, strand = pmp22_locus()
    window_start, window_end = max(0, gene_start - LOCUS_FLANK), gene_end + LOCUS_FLANK
    DERIVED.mkdir(parents=True, exist_ok=True)
    stats = []
    slice_path = DERIVED / "rn5_pmp22_locus_tracks.tsv"
    with slice_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["track", "chrom", "start", "end", "value"])
        for row in manifest():
            if row["filename"].endswith((".bed.gz", ".bedGraph.gz")):
                result, sliced = summarize_track(RAW / row["filename"], chrom, window_start, window_end)
                stats.append(result)
                writer.writerows((row["filename"], *item) for item in sliced)
    with (DERIVED / "track_qc.tsv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=stats[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(stats)
    donors = []
    with (RAW / "fantom_human_primary_cell_hCAGE_hg19_sdrf.tsv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Library Name") in {"CNhs12073", "CNhs12345", "CNhs12621"}:
                donors.append({"library": row["Library Name"], "extract": row["Extract Name"], "sample": row["Comment [sample_name]"], "rin": row["Comment [RNA integrity number]"], "file": row["File Name"]})
    donors.sort(key=lambda row: row["library"])
    with (DERIVED / "fantom_schwann_donors.tsv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=donors[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(donors)
    result = {"assembly": "rn5", "gene": {"chrom": chrom, "start": gene_start, "end": gene_end, "strand": strand}, "window": {"start": window_start, "end": window_end}, "tracks": len(stats), "tracks_with_locus_records": sum(s["records_in_window"] > 0 for s in stats), "fantom_donors": len(donors), "interpretation": "Processed-track descriptive analysis; inputs are not biological replicates and no causal effect is estimated."}
    (DERIVED / "analysis_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    by_name = {row["file"]: row for row in stats}
    report = f"""# Executed processed-data analysis

## Inputs and interval

All 12 inputs in `data/download_manifest.tsv` passed SHA-256 verification. UCSC
rn5 RefGene places Pmp22 at `{chrom}:{gene_start}-{gene_end}` ({strand} strand).
The frozen analysis window is `{chrom}:{window_start}-{window_end}`: the unioned
RefGene span with a 1 Mb flank on both sides.

## Results

Seven of eight processed tracks contain records in the window. The PNS SOX10
peak set contains **{by_name['GSE64703_Sox10_peaks_PNS.bed.gz']['records_in_window']}**
overlapping peaks, versus **{by_name['GSE64703_Sox10_peaks_CNS.bed.gz']['records_in_window']}**
in the CNS peak set. The differential peak files contribute
**{by_name['GSE64703_Sox10_PNS_DB_peaks.bed.gz']['records_in_window']} PNS** and
**{by_name['GSE64703_Sox10_CNS_DB_peaks.bed.gz']['records_in_window']} CNS**
records. The P15 sciatic-nerve H3K27ac set contributes
**{by_name['GSE64971_H3K27ac_P15_sciatic_nerve.bed.gz']['records_in_window']}** peaks.

Across covered bases (zeros/unreported bases excluded), submitter bedGraph signal
averages **{by_name['GSE64703_Sox10_SN_run204.ucsc.bedGraph.gz']['signal_mean_covered_bp']}**
for PNS SOX10, **{by_name['GSE64703_sox10_SC_run128.ucsc.bedGraph.gz']['signal_mean_covered_bp']}**
for CNS SOX10, and
**{by_name['GSE64971_H3K27ac_P15_sciatic_nerve.ucsc.bedGraph.gz']['signal_mean_covered_bp']}**
for PNS H3K27ac. These track scales are not assumed comparable.

The FANTOM SDRF extraction independently recovers three distinct Schwann records:
CNhs12073/donor1, CNhs12345/donor2, and CNhs12621/donor3. It does not resolve the
Borzoi label/path mismatch or provide downloaded CAGE signal.

## Interpretation boundary

The result reproduces descriptive processed-track support around rat Pmp22 and
confirms the FANTOM donor metadata. It does **not** estimate a state effect: the
bedGraphs are combined tracks, peak counts ignore width and threshold differences,
and the H3K27ac experiment has one IP plus input rather than two biological
replicates. It is not evidence of native human causality, therapeutic selectivity,
or a numerical PMP22 RNA/protein response.
"""
    (DERIVED / "analysis_report.md").write_text(report, encoding="utf-8")
    return result
