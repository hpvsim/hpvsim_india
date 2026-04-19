"""
Side-by-side comparison of fig2 plot-ready baselines across hpvsim versions.

Expects each baseline dir under `results/` to contain:
  - fig2_age_causal_summary.csv
  - fig2_dwelltime_summary.csv
  - fig2_outcomes.csv

Usage:
  python compare_fig2.py --baselines v2.2.6_baseline v2.3.0_baseline
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _bxp_stats(row, label):
    return dict(med=row['med'], q1=row['q1'], q3=row['q3'],
                whislo=row['whislo'], whishi=row['whishi'],
                fliers=[], label=label)


def compare_panel_ab(ax, dfs_by_version, label_col, title):
    """Grouped boxplots from summary CSVs."""
    labels = list(dfs_by_version[next(iter(dfs_by_version))][label_col])
    n_versions = len(dfs_by_version)
    n_labels = len(labels)
    group_width = 0.8
    box_width = group_width / n_versions

    colors = plt.cm.tab10(np.linspace(0, 1, n_versions))
    positions = []
    all_stats = []
    for vi, (version, df) in enumerate(dfs_by_version.items()):
        offset = (vi - (n_versions - 1) / 2) * box_width
        for li, label in enumerate(labels):
            row = df[df[label_col] == label].iloc[0]
            all_stats.append(_bxp_stats(row, f'{label}\n({version})' if vi == 0 else ''))
            positions.append(li + offset)

    bp = ax.bxp(all_stats, positions=positions, widths=box_width * 0.9,
                patch_artist=True, showfliers=False, manage_ticks=False)
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i // n_labels])
        patch.set_alpha(0.7)

    ax.set_xticks(np.arange(n_labels))
    ax.set_xticklabels(labels)
    ax.set_title(title)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=colors[vi], alpha=0.7, label=v)
                       for vi, v in enumerate(dfs_by_version)],
              loc='best', fontsize=9)


def compare_panel_cd(ax, outcomes_by_version, layers, end_years, title, denom_kind):
    """Overlay line plots of stacked-region proportions across versions."""
    linestyles = ['-', '--', ':', '-.']
    for vi, (version, outcomes) in enumerate(outcomes_by_version.items()):
        total_alive = outcomes['total']
        if denom_kind == 'all':
            layer_data = {
                'Cleared': outcomes['cleared'] / total_alive * 100,
                'Infection w/o HSIL': outcomes['persisted'] / total_alive * 100,
                'HSIL': (outcomes['progressed'] + outcomes['cancer'] + outcomes['dead']) / total_alive * 100,
            }
        else:
            total_persisted = total_alive - outcomes['cleared']
            layer_data = {
                'Infection w/o HSIL': outcomes['persisted'] / total_persisted * 100,
                'HSIL': outcomes['progressed'] / total_persisted * 100,
                'Cancer': (outcomes['cancer'] + outcomes['dead']) / total_persisted * 100,
            }

        dt_years = outcomes['years'][1] - outcomes['years'][0]
        end_ind = int(1 / dt_years) * end_years
        xs = outcomes['years'][:end_ind]
        ls = linestyles[vi % len(linestyles)]
        colors = ['#2db5a7', '#eddc42', '#e67f2c', '#871a6c']
        for li, label in enumerate(layers):
            offset_color = 0 if denom_kind == 'all' else 1
            ax.plot(xs, layer_data[label][:end_ind], color=colors[li + offset_color],
                    linestyle=ls, label=f'{label} ({version})')

    ax.set_title(title)
    ax.set_xlabel('Time since infection')
    ax.set_ylabel('%')
    ax.legend(loc='best', fontsize=8, ncol=2)


def load_baseline(resfolder, baseline):
    d = Path(resfolder) / baseline
    try:
        return dict(
            ac=pd.read_csv(d / 'fig2_age_causal_summary.csv'),
            dw=pd.read_csv(d / 'fig2_dwelltime_summary.csv'),
            outcomes=pd.read_csv(d / 'fig2_outcomes.csv'),
        )
    except FileNotFoundError as e:
        print(f'Skipping {baseline}: {e}')
        return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baselines', nargs='+', default=['v2.2.6_baseline'])
    parser.add_argument('--resfolder', default='results')
    parser.add_argument('--outpath', default='figures/compare/fig2_compare.png')
    args = parser.parse_args()

    data = {}
    for b in args.baselines:
        loaded = load_baseline(args.resfolder, b)
        if loaded is not None:
            data[b] = loaded
    if not data:
        raise SystemExit('No baselines with fig2 data found.')

    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), layout='tight')

    compare_panel_ab(axes[0, 0], {v: d['ac'] for v, d in data.items()},
                     label_col='Health event',
                     title='(A) Age distribution of key health events')
    axes[0, 0].set_ylabel('Age')
    axes[0, 0].set_ylim([0, 100])

    compare_panel_ab(axes[0, 1], {v: d['dw'] for v, d in data.items()},
                     label_col='state', title='(B) Dwelltimes')
    axes[0, 1].set_ylabel('Dwelltime')

    compare_panel_cd(axes[1, 0], {v: d['outcomes'] for v, d in data.items()},
                     layers=['Cleared', 'Infection w/o HSIL', 'HSIL'],
                     end_years=10, title='(C) All outcomes', denom_kind='all')

    compare_panel_cd(axes[1, 1], {v: d['outcomes'] for v, d in data.items()},
                     layers=['Infection w/o HSIL', 'HSIL', 'Cancer'],
                     end_years=40, title='(D) Conditional on persistence', denom_kind='persisted')

    fig.savefig(args.outpath, dpi=120)
    print(f'Wrote {args.outpath}')
