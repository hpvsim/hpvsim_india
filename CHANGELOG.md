# Changelog

## 2026-04-18 — HPVsim v2.2.6 lift

- Split every plot script into a VM-side `--run-sim` step that saves lightweight CSVs and a local plot step that reads them.
- Added `compare_fig2.py` for cross-version comparison of Fig 2.
- Froze plot-ready baselines under `results/v2.2.6_baseline/` (all figures) and `results/v1.0.0_published/` (figS1 + figS2 panels B/D/E, extracted from original paper `.obj` artifacts).
- Dropped MySQL `storage` requirement from calibration; in-memory Optuna by default.
- Fixed figS1 Panel B: boxplots and DHS scatter now share the AgeRange axis.
- Added `.gitignore` excluding raw per-event extractions, `india_calib.obj`, and versioned `.obj`/`.msim` artifacts.
- Refreshed README with the new workflow and pinned HPVsim version.

## 2024-01 — v2.0 era

- Original India calibration and sexual-behavior analysis for the HPVsim methods manuscript (Fig 2 and Figs S1–S3).
