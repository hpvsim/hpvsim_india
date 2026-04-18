"""
Plot calibration to India.

Two modes:
  python plot_figS2.py --run-sim   # run calibration + extract CSVs (VM, slow)
  python plot_figS2.py             # plot from saved CSVs (local)
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc
import seaborn as sns

import run_sim as rs
import utils as ut

do_show = True


# %% Data extraction (run on VM after calibration)
def save_figS2_data(calib, res_to_plot=100, resfolder='results'):
    indices = calib.df.iloc[0:res_to_plot, 0].values.astype(int)

    rows = []
    for i in indices:
        cancers = calib.analyzer_results[i]['cancers'][2020]
        for bin_idx, val in enumerate(cancers):
            rows.append({'trial_idx': int(i), 'bin': bin_idx, 'value': float(val)})
    pd.DataFrame(rows).to_csv(f'{resfolder}/figS2_cancers_by_age.csv', index=False)

    for rkey in ['cin_genotype_dist', 'cancerous_genotype_dist']:
        rows = []
        for i in indices:
            run = calib.sim_results[i][rkey]
            arr = [run] if sc.isnumber(run) else list(run)
            for bin_idx, val in enumerate(arr):
                rows.append({'trial_idx': int(i), 'bin': bin_idx, 'value': float(val)})
        pd.DataFrame(rows).to_csv(f'{resfolder}/figS2_{rkey}.csv', index=False)

    calib.target_data[0].to_csv(f'{resfolder}/figS2_target_cancers.csv', index=False)
    calib.target_data[1].to_csv(f'{resfolder}/figS2_target_cin_genotype.csv', index=False)
    calib.target_data[2].to_csv(f'{resfolder}/figS2_target_cancerous_genotype.csv', index=False)


# %% Plotting (run locally)
def plot_calib(resfolder='results', outpath='figures/figS2.png'):
    ut.set_font(size=15)
    fig = plt.figure(layout='tight', figsize=(10, 7))
    prev_col = '#5f5cd2'
    canc_col = '#c1981d'
    ms = 80
    gen_cols = sc.gridcolors(4)

    gs0 = fig.add_gridspec(2, 1)
    gs00 = gs0[0].subgridspec(1, 1)
    gs01 = gs0[1].subgridspec(1, 4)

    # ###############
    # # Panel A: HPV prevalence by age
    # ###############
    res = sc.loadobj(f'{resfolder}/india_msim.obj')
    year = 2020
    ind = sc.findinds(res['year'], year)[0]

    pre_cins = {}
    ts = 0.5
    for which in ['values', 'low', 'high']:
        this_res = res['n_precin_by_age'][which][:, ind]
        pre_cins[which] = [
            sum(this_res[3:5]) / sum(res['n_females_alive_by_age'][3:5, ind]),
            sum(this_res[5:7]) * ts / sum(res['n_females_alive_by_age'][5:7, ind]),
            sum(this_res[7:9]) * ts / sum(res['n_females_alive_by_age'][7:9, ind]),
            sum(this_res[9:11]) * ts / sum(res['n_females_alive_by_age'][9:11, ind]),
            sum(this_res[11:13]) * ts / sum(res['n_females_alive_by_age'][11:13, ind]),
        ]

    ax = fig.add_subplot(gs01[:2])
    age_labels = ['15-25', '25-34', '35-44', '45-54', '55-64']
    x = np.arange(len(age_labels))
    best = np.array([.15, .13, .13, .13, .13])
    lowererr = np.array([0.025, 0.015, 0.02, 0.025, 0.08])
    uppererr = np.array([0.02, 0.01, 0.015, 0.03, 0.09])
    err = [lowererr, uppererr]

    ax.plot(x, pre_cins['values'], color=prev_col)
    ax.fill_between(x, pre_cins['low'], pre_cins['high'], color=prev_col, alpha=0.3)
    ax.errorbar(x, best, yerr=err, ls='none', marker='d', markersize=ms / 10, color='k')
    ax.set_ylim([0, 0.25])
    ax.set_xticks(x, age_labels)
    ax.set_xlabel('Age')
    ax.set_title('Detectable HPV prevalence,\n normal cervical cytology, 2020')

    ###############
    # Panel B: Cancers by age (from calib extract)
    ###############
    ax = fig.add_subplot(gs00[0])
    cancers_df = pd.read_csv(f'{resfolder}/figS2_cancers_by_age.csv')
    target_cancers = pd.read_csv(f'{resfolder}/figS2_target_cancers.csv')
    age_labels = ['0', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85']
    x = np.arange(len(age_labels))

    sns.lineplot(ax=ax, x='bin', y='value', data=cancers_df, color=canc_col, errorbar=('pi', 95))
    ax.scatter(x, target_cancers['value'].values, marker='d', s=ms, color='k')
    ax.set_ylim([0, 20_000])
    ax.set_xticks(x, age_labels)
    ax.set_xlabel('Age')
    ax.set_title('Cancers by age, 2020')

    ###############
    # Panels D, E: CIN + cancer by genotype
    ###############
    for ai, rkey in enumerate(['cin_genotype_dist', 'cancerous_genotype_dist']):
        ax = fig.add_subplot(gs01[ai + 2])
        model_df = pd.read_csv(f'{resfolder}/figS2_{rkey}.csv')
        target_name = {'cin_genotype_dist': 'cin_genotype', 'cancerous_genotype_dist': 'cancerous_genotype'}[rkey]
        target_df = pd.read_csv(f'{resfolder}/figS2_target_{target_name}.csv')

        sns.boxplot(ax=ax, x='bin', y='value', data=model_df, palette=gen_cols, showfliers=False)
        ax.scatter(np.arange(len(target_df)), target_df['value'].values, color='k', marker='d', s=ms)
        ax.set_ylim([0, 1])
        ax.set_xticks(np.arange(4), ['16', '18', 'Hi5', 'OHR'])
        ax.set_title({'cin_genotype_dist': 'CIN2+ by genotype',
                      'cancerous_genotype_dist': 'Cancers by genotype'}[rkey])

    fig.tight_layout()
    sc.savefig(outpath, dpi=300)
    if do_show:
        plt.show()


# %% Run as a script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run calibration and save CSVs (heavy, VM-side)')
    parser.add_argument('--resfolder', default='results')
    parser.add_argument('--figpath', default='figures/figS2.png')
    parser.add_argument('--res-to-plot', type=int, default=100)
    args = parser.parse_args()

    if args.run_sim:
        sim, calib = rs.run_calib(n_trials=rs.n_trials, n_workers=rs.n_workers,
                                  do_save=True, filestem='')
        save_figS2_data(calib, res_to_plot=args.res_to_plot, resfolder=args.resfolder)
        print(f'Saved figS2 CSVs to {args.resfolder}')
    else:
        plot_calib(resfolder=args.resfolder, outpath=args.figpath)
        print('Done.')
