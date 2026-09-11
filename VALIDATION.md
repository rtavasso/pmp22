# Release validation

- 28 unit/integrity tests pass in the analysis checkout and a clean Git clone.
- All 112 locked public resources and the original 12-input manifest verify.
- A fresh checkout restores all 62 public metadata snapshots with network calls
  disabled. These preserve exact API bytes despite changing remote metadata.
- The review correction reruns human mapping/activity, bulk context/control
  selection, effect curation, dossiers, matrix compilation and report generation.
  The comparability analysis recomputes CpG from verified reference sequences
  and scores frozen predictions in explicit TSS strata.
- The original split/protocol remain byte-identical. Predictions, coefficients
  and primary metrics retain their original Git content (working-tree text
  line endings may differ). Scores recompute with round-trip float parsing.
  The original full benchmark replay is retained as prior validation history;
  this correction does not refit or reselect the original models.
- The revised two-panel scientific benchmark figure was visually inspected.
  The standalone HTML was regenerated and its anchors, embedded tables and
  downloadable payloads checked. The original browser/figure inspection is
  retained in prior validation history.
- Git whitespace checks pass with the frozen CRLF split explicitly preserved.

The tested code commit and timestamp are recorded in data/validation.json.
These checks establish software/data integrity, not native human causality,
human-state generalization, clinical efficacy or publication novelty.
