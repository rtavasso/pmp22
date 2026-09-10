# Benchmark contract

No benchmark score is fabricated in this seed release.  A run is admissible only
when `splits/`, exact labels, source-independent units, model/checkpoint versions,
input/output transforms, seeds, and the exposure audit are frozen.

Required tasks are kept separate:

* sequence/accessibility → held-out accessibility profile;
* sequence → held-out RNA/CAGE profile; and
* sequence intervention → an experimentally measured intervention effect.

The first two do not imply the third.  Each task compares a distance-only
baseline, motifs, conservation, and measured accessibility when that information
is available to the selected model on identical examples.  Primary reporting
includes confidence intervals, calibration, source-level resampling, robustness
to region definitions, and null findings.

Stop rather than score when the test interval occurred in the selected training
fold and no independent endpoint exists, when the apparent state split is a
study split, or when labels lack sufficient independent elements.  A future run
must write immutable split IDs under `benchmark/splits/` and outputs under an
ignored run directory.

