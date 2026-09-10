import csv
import json
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path

from pmp22_atlas.cli import DATA, SCHEMAS, build, read_table, validate
from pmp22_atlas.acquire import manifest


class AtlasTest(unittest.TestCase):
    def test_html_report_is_standalone_and_contains_executed_results(self):
        report = (DATA.parent / "report.html").read_text(encoding="utf-8")
        parser = HTMLParser()
        parser.feed(report)
        self.assertIn("<strong>35</strong>", report)
        self.assertIn("<strong>3</strong>", report)
        self.assertIn("GSE64703", report)
        self.assertNotIn("src=\"http", report)

    def test_all_tables_validate_and_have_rows(self):
        counts = validate()
        self.assertEqual(set(counts), set(SCHEMAS))
        self.assertTrue(all(count > 0 for count in counts.values()))

    def test_unknown_coordinates_are_not_fabricated(self):
        for row in read_table("regions.tsv"):
            values = (row["chrom"], row["start"], row["end"])
            self.assertIn(sum(bool(value) for value in values), (0, 3))

    def test_fantom_donors_are_distinct(self):
        rows = [r for r in read_table("samples.tsv") if r["source_id"] == "FANTOM5"]
        self.assertEqual(3, len(rows))
        self.assertEqual(3, len({r["sample_id"] for r in rows}))
        self.assertEqual(3, len({r["donor_or_pool"] for r in rows}))

    def test_build_is_complete_and_deterministic(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            summary = build(Path(first))
            build(Path(second))
            expected = {"evidence_matrix.tsv", "candidate_panel.tsv", "source_overlap.dot", "atlas_summary.json"}
            self.assertEqual(expected, {p.name for p in Path(first).iterdir()})
            for name in expected:
                self.assertEqual((Path(first) / name).read_bytes(), (Path(second) / name).read_bytes())
            self.assertEqual(5, summary["candidate_count"])
            self.assertEqual("provisional", json.loads((Path(first) / "atlas_summary.json").read_text())["status"])

    def test_effect_region_foreign_keys(self):
        regions = {r["region_id"] for r in read_table("regions.tsv")}
        self.assertTrue(all(r["region_id"] in regions for r in read_table("assay_effects.tsv")))

    def test_download_manifest_is_pinned(self):
        rows = manifest()
        self.assertEqual(12, len(rows))
        self.assertTrue(all(row["url"].startswith("https://") for row in rows))
        self.assertTrue(all(len(row["sha256"]) == 64 for row in rows))

    def test_committed_analysis_has_expected_locus_results(self):
        summary = json.loads((DATA / "derived" / "analysis_summary.json").read_text())
        self.assertEqual({"chrom": "chr10", "start": 49316999, "end": 49346994, "strand": "+"}, summary["gene"])
        self.assertEqual(7, summary["tracks_with_locus_records"])
        with (DATA / "derived" / "track_qc.tsv").open(newline="") as handle:
            rows = {row["file"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual("35", rows["GSE64703_Sox10_peaks_PNS.bed.gz"]["records_in_window"])
        self.assertEqual("3", rows["GSE64703_Sox10_peaks_CNS.bed.gz"]["records_in_window"])


if __name__ == "__main__":
    unittest.main()
