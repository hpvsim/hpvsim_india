# HPVsim model for India

Code for creating and calibrating a model of HPV transmission and progression in India, and for reproducing the India-specific figures of the HPVsim methods manuscript.

**Results in this repository were produced with HPVsim v2.2.6.** Baselines for comparison across hpvsim versions live in [`results/v2.2.6_baseline/`](results/v2.2.6_baseline/) (current) and [`results/v1.0.0_published/`](results/v1.0.0_published/) (extracted from the original paper artifacts, where available).

## Installation

```bash
pip install hpvsim==2.2.6 seaborn
```

Python 3.9+.

## Workflow: heavy sims on VM, plots locally

Every plot script has two modes: `--run-sim` runs the heavy calculation and saves lightweight CSVs under `results/`; running without flags loads the CSVs and produces the figure. The intended flow is:

1. **VM:** run `python plot_XXX.py --run-sim`, commit and push the resulting `results/*.csv` files
2. **Local:** pull, run `python plot_XXX.py` to render the figure from the CSVs

Per-event extractions (e.g. age at cancer onset for a million agents) are not committed — only boxplot summary statistics are. See `.gitignore` for the exclusion list.

### Plotting scripts

| Script | Figure | Heavy step | Key artifacts |
|---|---|---|---|
| `plot_fig2.py` | Fig 2 — implied natural history | 1M-agent sim with custom analyzers | `fig2_age_causal_summary.csv`, `fig2_dwelltime_summary.csv`, `fig2_outcomes.csv` |
| `plot_figS1.py` | Fig S1 — sexual behavior | Short sim with AFS / prop-married / snapshot analyzers | `model_sb_AFS.csv`, `model_sb_prop_married.csv`, `model_age_diffs.csv`, `model_casual.csv`, `partners.obj` |
| `plot_figS2.py` | Fig S2 — calibration | Full calibration (10k trials) | `figS2_cancers_by_age.csv`, `figS2_cin_genotype_dist.csv`, `figS2_cancerous_genotype_dist.csv` + 3 target CSVs |
| `plot_figS3.py` | Fig S3 — India cancer burden vs Globocan | Pre-computed `india_msim.obj` | `india_msim.obj` |

### Other scripts

- `run_sim.py`: sim-building, calibration, sexual-behavior extraction (called by the `--run-sim` modes)
- `save_figS3_baseline.py`: freezes figS3 plot-ready CSVs into a versioned baseline dir
- `compare_fig2.py`: generates side-by-side fig2 comparison across baseline dirs
- `utils.py`, `behavior_inputs.py`: supporting utilities and parameters

## Cross-version comparison

Each baseline lives under `results/<version>/` and contains CSVs (never pickles). To produce a side-by-side fig2 comparison:

```bash
python compare_fig2.py --baselines v2.2.6_baseline v2.3.0_baseline
```

### Adding a future-version baseline (v2.3, v3.0, ...)

When a new HPVsim version ships:

```bash
# 1. On VM, in a clean env pinned to the new version
conda create -n hpvsim230 python=3.11 -y && conda activate hpvsim230
pip install hpvsim==2.3.0 seaborn optuna
# 2. Re-run each script that has a --run-sim mode
python plot_fig2.py --run-sim       # heavy: 1M-agent sim
python plot_figS1.py --run-sim      # sexual-behavior extraction
python plot_figS2.py --run-sim      # full calibration (slowest)
# (plot_figS3 loads the committed india_msim.obj; regenerate via run_sim.py's run_parsets flow)
# 3. Freeze the fresh CSVs into a versioned baseline dir
mkdir -p results/v2.3.0_baseline
cp results/*.csv results/v2.3.0_baseline/
# 4. Commit + push + compare
python compare_fig2.py --baselines v2.2.6_baseline v2.3.0_baseline
```

The comparison script extends trivially to v3.0 by appending the new baseline name to `--baselines`.

## Inputs

- `data/` — input data files (DHS debut/marriage, Globocan, IARC targets)
- `behavior_inputs.py` — additional behavioral parameters

## Further information

See [hpvsim.org](https://hpvsim.org) and [docs.hpvsim.org](https://docs.hpvsim.org).
