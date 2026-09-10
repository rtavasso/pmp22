"""Validate and compile the versioned Project A manifests."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

SCHEMAS = {
    "regions.tsv": ("region_id", "assembly", "chrom", "start", "end", "strand", "element_class", "source_id", "evidence_status", "orthology_status"),
    "samples.tsv": ("sample_id", "source_id", "species", "tissue", "state", "assay", "role", "donor_or_pool", "download_status", "permission"),
    "source_overlap.tsv": ("parent_id", "child_id", "relationship", "independence", "notes"),
    "assay_effects.tsv": ("effect_id", "region_id", "source_id", "comparison", "outcome", "effect_value", "effect_unit", "replicates", "evidence_layer", "value_origin", "species", "cell_context", "tested_interval", "control", "uncertainty", "locator", "limitations"),
    "claim_ledger.tsv": ("claim_id", "claim", "unit", "estimand", "source_ids", "independence", "limiting_assumptions", "allowed_conclusion"),
    "model_exposure.tsv": ("model_id", "checkpoint", "source_id", "interval_exposure", "source_exposure", "assay_exposure", "status", "notes"),
    "software_lock.tsv": ("component", "version", "reference", "purpose"),
    "prior_study_inventory.tsv": ("source_id", "citation", "species", "assay", "result_scope", "access_status", "novelty_note"),
}


def read_table(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SCHEMAS[name]:
            raise ValueError(f"{name}: expected columns {SCHEMAS[name]}, got {reader.fieldnames}")
        return list(reader)


def validate() -> dict[str, int]:
    tables = {name: read_table(name) for name in SCHEMAS}
    regions = tables["regions.tsv"]
    ids = [r["region_id"] for r in regions]
    if len(ids) != len(set(ids)):
        raise ValueError("regions.tsv: duplicate region_id")
    for row in regions:
        coords = (row["chrom"], row["start"], row["end"])
        if any(coords) and not all(coords):
            raise ValueError(f"{row['region_id']}: coordinates must be wholly known or empty")
        if all(coords):
            try:
                start, end = int(row["start"]), int(row["end"])
            except ValueError as exc:
                raise ValueError(f"{row['region_id']}: non-integer coordinate") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{row['region_id']}: invalid half-open interval")
        if row["evidence_status"] not in {"measured", "predicted", "proposed", "missing"}:
            raise ValueError(f"{row['region_id']}: invalid evidence_status")
    region_ids = set(ids)
    for effect in tables["assay_effects.tsv"]:
        if effect["region_id"] not in region_ids:
            raise ValueError(f"{effect['effect_id']}: unknown region {effect['region_id']}")
    return {name: len(rows) for name, rows in tables.items()}


def write_tsv(path: Path, rows: list[dict[str, str]], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(output: Path) -> dict[str, object]:
    counts = validate()
    output.mkdir(parents=True, exist_ok=True)
    regions = read_table("regions.tsv")
    effects = read_table("assay_effects.tsv")
    by_region: dict[str, list[dict[str, str]]] = {}
    for effect in effects:
        by_region.setdefault(effect["region_id"], []).append(effect)
    matrix = []
    for region in regions:
        linked = by_region.get(region["region_id"], []) or [{}]
        for effect in linked:
            matrix.append({
                "region_id": region["region_id"], "element_class": region["element_class"],
                "assembly": region["assembly"], "interval": (f"{region['chrom']}:{region['start']}-{region['end']}" if region["chrom"] else "UNRESOLVED"),
                "orthology_status": region["orthology_status"], "evidence_status": region["evidence_status"],
                "effect_id": effect.get("effect_id", ""), "outcome": effect.get("outcome", ""),
                "evidence_layer": effect.get("evidence_layer", "missing"),
                "source_id": effect.get("source_id",region["source_id"]),
                "species": effect.get("species","human reference annotation"),
                "value": effect.get("effect_value",""),"unit": effect.get("effect_unit",""),
                "tested_interval": effect.get("tested_interval","annotated reference interval"),
                "limitation": effect.get("limitations","Annotation is not a measured regulatory effect"),
            })
    region_map={r["region_id"]:r for r in regions}
    def append_observation(rid,source,layer,value,unit,species,limitation):
        r=region_map[rid]
        matrix.append(dict(region_id=rid,element_class=r["element_class"],assembly=r["assembly"],
            interval=f'{r["chrom"]}:{r["start"]}-{r["end"]}',orthology_status=r["orthology_status"],
            evidence_status="measured" if species.startswith("human") else "proposed",
            effect_id="",outcome="interval overlap or promoter initiation",evidence_layer=layer,
            source_id=source,species=species,value=str(value),unit=unit,tested_interval="observational reference overlap",
            limitation=limitation))
    def optional_table(name):
        path=DATA/"human"/name
        if not path.exists():return []
        with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
    for r in optional_table("region_activity.tsv"):
        for field,source in [("schwann_atac_unique_peak_intervals","ENCFF602YVR"),("broad_atac_unique_peak_intervals","ENCFF632FMT")]:
            append_observation(r["region_id"],source,"accessibility",r[field],"unique peak intervals","human pooled Schwann",
                "Threshold-specific observation from same two-donor pool; not a causal expression effect")
        append_observation(r["region_id"],"UCSC phastCons100way","conservation",r["phastcons100way_mean"],
            "mean phastCons over covered bases","human reference / multispecies alignment","Conservation does not establish enhancer activity")
    for r in optional_table("promoter_cage.tsv"):
        if r["upstream_flank_bp"]!="100":continue
        append_observation(r["promoter"],r["donor"],"CAGE",r["minus_strand_tags"],"minus-strand tags",
            "human cultured Schwann","First-exon-family window plus 100 bp flank; not mature-RNA isoform abundance")
    for r in optional_table("region_bulk_support.tsv"):
        append_observation(r["region_id"],r["file"],r["target"],r["unique_peak_intervals"],"unique peak intervals",
            "human bulk tibial nerve","Shared donors across assays; cell composition unresolved")
    fields = tuple(matrix[0])
    write_tsv(output / "evidence_matrix.tsv", matrix, fields)

    candidates = []
    dossier_dir = ROOT / "candidate_dossiers"
    for path in sorted(dossier_dir.glob("*.md")):
        meta = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith("- **") and ":** " in line:
                key, value = line[4:].split(":** ", 1)
                meta[key.lower()] = value
        candidates.append({"candidate": path.stem, "priority": meta["priority"], "direction": meta["direction"], "decisive_experiment": meta["decisive experiment"], "dossier": str(path.relative_to(ROOT))})
    candidates.sort(key=lambda row: (int(row["priority"]), row["candidate"]))
    write_tsv(output / "candidate_panel.tsv", candidates, tuple(candidates[0]))

    overlap = read_table("source_overlap.tsv")
    dot = ["digraph source_overlap {", '  rankdir="LR";']
    for row in overlap:
        label = row["relationship"].replace('"', "'")
        style = "dashed" if row["independence"] != "independent" else "solid"
        dot.append(f'  "{row["parent_id"]}" -> "{row["child_id"]}" [label="{label}", style="{style}"];')
    dot.append("}")
    (output / "source_overlap.dot").write_text("\n".join(dot) + "\n", encoding="utf-8")
    summary = {"status": "public_data_release; native_human_validation_pending", "table_rows": counts, "regions_with_coordinates": sum(bool(r["chrom"]) for r in regions), "regions_unresolved": sum(not bool(r["chrom"]) for r in regions), "candidate_count": len(candidates), "claim_boundary": "No new causal, clinical, or native-human perturbation result is asserted."}
    (output / "atlas_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pmp22_atlas")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate", help="validate all committed manifests")
    build_parser = sub.add_parser("build", help="compile deterministic atlas outputs")
    build_parser.add_argument("--output", type=Path, default=ROOT / "build")
    download_parser = sub.add_parser("download", help="download and checksum pinned public inputs")
    download_parser.add_argument("--force", action="store_true")
    sub.add_parser("analyze", help="analyze downloaded processed tracks at the rn5 Pmp22 locus")
    for command in ['human','context','benchmark','supplemental','report']:
        sub.add_parser(command, help=f'run {command} analysis (analysis dependencies required)')
    args = parser.parse_args(argv)
    if args.command in {"download", "analyze"}:
        from .acquire import analyze, download
        result = download(args.force) if args.command == "download" else analyze()
    elif args.command in {'human','context','benchmark','supplemental','report'}:
        import importlib
        importlib.import_module('pmp22_atlas.'+args.command).main()
        return 0
    else:
        result = validate() if args.command == "validate" else build(args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
