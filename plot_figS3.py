"""
Plot India cancer burden: HPVsim vs Globocan.

Loads plot-ready CSVs from a baseline dir. To regenerate from a fresh msim
run, use save_figS3_baseline.py --outdir results/v<version>_baseline.
"""
import argparse

import matplotlib.pyplot as plt
import pandas as pd
import sciris as sc

import utils as ut

do_show = True


def plot_india_burden(resfolder, outpath='figures/figS3.png'):
    ut.set_font(size=16)
    hpvsim_df = pd.read_csv(f'{resfolder}/figS3_hpvsim.csv')
    globocan_df = pd.read_csv(f'{resfolder}/figS3_globocan.csv')
    colors = sc.gridcolors(2)

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.plot(hpvsim_df['year'], hpvsim_df['asr_cancer_incidence'],
            color=colors[0], label='HPVsim', marker='o')
    ax.fill_between(hpvsim_df['year'], hpvsim_df['asr_cancer_incidence_low'],
                    hpvsim_df['asr_cancer_incidence_high'],
                    color=colors[0], alpha=0.3)
    ax.plot(globocan_df['year'], globocan_df['asr_cancer_incidence'],
            marker='s', color=colors[1], label='Globocan')
    ax.set_ylabel('Age-standardized cancer incidence (per 100k)')
    ax.legend()
    ax.set_ylim([0, 25])
    fig.tight_layout()
    sc.savefig(outpath, dpi=100)
    if do_show:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/figS3.png')
    args = parser.parse_args()
    plot_india_burden(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
