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

import run_sim as rs
import utils as ut

do_show = True


# %% Data extraction (run on VM after calibration)
def save_figS2_data(calib, res_to_plot=100, resfolder='results'):
    indices = calib.df.iloc[0:res_to_plot, 0].values.astype(int)

    # Summarize cancers-by-age across trials: median + 95% percentile interval per bin
    per_trial = np.array([calib.analyzer_results[i]['cancers'][2020] for i in indices])
    pd.DataFrame({
        'bin': np.arange(per_trial.shape[1]),
        'median': np.median(per_trial, axis=0),
        'pi95_low': np.percentile(per_trial, 2.5, axis=0),
        'pi95_high': np.percentile(per_trial, 97.5, axis=0),
    }).to_csv(f'{resfolder}/figS2_cancers_by_age.csv', index=False)

    # Summarize genotype distributions to boxplot stats per bin
    for rkey in ['cin_genotype_dist', 'cancerous_genotype_dist']:
        per_trial = []
        for i in indices:
            run = calib.sim_results[i][rkey]
            per_trial.append([run] if sc.isnumber(run) else list(run))
        per_trial = np.array(per_trial)
        rows = []
        for bi in range(per_trial.shape[1]):
            arr = per_trial[:, bi]
            q1, med, q3 = np.percentile(arr, [25, 50, 75])
            iqr = q3 - q1
            lo = float(arr[arr >= q1 - 1.5 * iqr].min())
            hi = float(arr[arr <= q3 + 1.5 * iqr].max())
            rows.append({'bin': bi, 'q1': float(q1), 'med': float(med), 'q3': float(q3),
                         'whislo': lo, 'whishi': hi})
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

    # Panel A: HPV prevalence by age, 2020
    prev_df = pd.read_csv(f'{resfolder}/figS2_hpv_prevalence.csv')
    x = np.arange(len(prev_df))
    best = np.array([.15, .13, .13, .13, .13])
    err = [np.array([0.025, 0.015, 0.02, 0.025, 0.08]),
           np.array([0.02, 0.01, 0.015, 0.03, 0.09])]

    ax = fig.add_subplot(gs01[:2])
    ax.plot(x, prev_df['values'], color=prev_col)
    ax.fill_between(x, prev_df['low'], prev_df['high'], color=prev_col, alpha=0.3)
    ax.errorbar(x, best, yerr=err, ls='none', marker='d', markersize=ms / 10, color='k')
    ax.set_ylim([0, 0.25])
    ax.set_xticks(x, prev_df['label'].tolist())
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

    cancers_df = cancers_df.sort_values('bin')
    ax.plot(cancers_df['bin'], cancers_df['median'], color=canc_col)
    ax.fill_between(cancers_df['bin'], cancers_df['pi95_low'],
                    cancers_df['pi95_high'], color=canc_col, alpha=0.3)
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

        model_df = model_df.sort_values('bin')
        bxp_stats = [dict(med=r.med, q1=r.q1, q3=r.q3,
                          whislo=r.whislo, whishi=r.whishi, fliers=[],
                          label=str(int(r.bin))) for _, r in model_df.iterrows()]
        bp = ax.bxp(bxp_stats, positions=model_df['bin'].values,
                    patch_artist=True, showfliers=False, manage_ticks=False)
        for patch, color in zip(bp['boxes'], gen_cols):
            patch.set_facecolor(color)
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
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline',
                        help='Dir with plot-ready CSVs (for plot mode only)')
    parser.add_argument('--figpath', default='figures/figS2.png')
    parser.add_argument('--res-to-plot', type=int, default=100)
    args = parser.parse_args()

    if args.run_sim:
        sim, calib = rs.run_calib(n_trials=rs.n_trials, n_workers=rs.n_workers,
                                  do_save=True, filestem='')
        save_figS2_data(calib, res_to_plot=args.res_to_plot, resfolder='results')
        print('Saved figS2 CSVs to results/ (copy to a versioned baseline dir to commit)')
    else:
        plot_calib(resfolder=args.resfolder, outpath=args.figpath)
        print('Done.')
