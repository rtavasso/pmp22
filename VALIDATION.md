# Release validation

- 24 unit/integrity tests pass in the analysis checkout and a clean Git clone.
- All 112 locked public resources and the original 12-input manifest verify.
- A fresh checkout restores all 62 public metadata snapshots with network calls
  disabled. These preserve exact API bytes despite changing remote metadata.
- The primary benchmark replay preserves its frozen split/protocol and
  reproduces the metrics. Prediction CSVs use round-trip float parsing to
  preserve tied scores when recalculating AP.
- Three scientific figures and the browser report were visually inspected;
  chart fonts have portable fallbacks.
- Git whitespace checks pass with the frozen CRLF split explicitly preserved.

The tested code commit and timestamp are recorded in data/validation.json.
These checks establish software/data integrity, not native human causality,
human-state generalization, clinical efficacy or publication novelty.
