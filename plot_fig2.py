"""
Plot implied natural history.

Two modes:
  python plot_fig2.py --run-sim   # run sim, save lightweight CSVs (VM)
  python plot_fig2.py             # plot from saved CSVs (local)
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc

import run_sim as rs
import utils as ut

do_show = True


# %% Data extraction (run on VM)
def _boxstats(values):
    arr = np.asarray(values)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    iqr = q3 - q1
    whislo = float(arr[arr >= q1 - 1.5 * iqr].min())
    whishi = float(arr[arr <= q3 + 1.5 * iqr].max())
    return dict(q1=float(q1), med=float(med), q3=float(q3),
                whislo=whislo, whishi=whishi, n=int(arr.size))


def get_age_causal_summary(sim):
    dt = sim.get_analyzer('age_causal')
    ages = np.array(dt.age_causal)
    categories = [
        ('Causal infection', np.array(dt.age_causal)[ages < 50]),
        ('CIN2+', np.array(dt.age_cin)[ages < 65]),
        ('Cancer', np.array(dt.age_cancer)[ages < 90]),
    ]
    rows = []
    for label, arr in categories:
        row = _boxstats(arr)
        row['Health event'] = label
        rows.append(row)
    return pd.DataFrame(rows)


def get_dwelltime_summary(sim):
    a = sim.get_analyzer('dwelltime_by_genotype')
    rows = []
    for cin in ['precin', 'cin']:
        dt = a.dwelltime[cin]
        arr = dt[0] + dt[1] + dt[2] + dt[3]
        row = _boxstats(arr)
        row['state'] = cin.upper()
        rows.append(row)
    return pd.DataFrame(rows)


def get_outcomes_df(sim):
    a = sim.get_analyzer('outcomes_by_year')
    res = a.results
    return pd.DataFrame({
        'years': a.durations,
        'total': res['total'],
        'cleared': res['cleared'],
        'persisted': res['persisted'],
        'progressed': res['progressed'],
        'cancer': res['cancer'],
        'dead': res['dead'],
    })


def save_fig2_data(sim, resfolder='results'):
    get_age_causal_summary(sim).to_csv(f'{resfolder}/fig2_age_causal_summary.csv', index=False)
    get_dwelltime_summary(sim).to_csv(f'{resfolder}/fig2_dwelltime_summary.csv', index=False)
    get_outcomes_df(sim).to_csv(f'{resfolder}/fig2_outcomes.csv', index=False)


# %% Plotting (run locally)
def _bxp_stats(df, label_col):
    return [
        dict(med=r['med'], q1=r['q1'], q3=r['q3'],
             whislo=r['whislo'], whishi=r['whishi'],
             fliers=[], label=r[label_col])
        for _, r in df.iterrows()
    ]


def plot_fig2(ac_summary, dw_summary, outcomes_df, outpath='figures/fig2.png'):
    ut.set_font(size=16)
    ac_colors = sc.gridcolors(3)
    stacked_colors = ['#2db5a7', '#eddc42', '#e67f2c', '#871a6c']
    dw_colors = stacked_colors[1:3]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()

    # (A) Age of causal events
    ax = axes[0]
    bp = ax.bxp(_bxp_stats(ac_summary, 'Health event'), patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], ac_colors):
        patch.set_facecolor(color)
    ax.set_title('(A) Age distribution of key health events')
    ax.set_xlabel('')
    ax.set_ylabel('Age')
    ax.set_ylim([0, 100])

    # (B) Dwelltimes
    ax = axes[1]
    bp = ax.bxp(_bxp_stats(dw_summary, 'state'), patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], dw_colors):
        patch.set_facecolor(color)
    ax.set_xlabel('')
    ax.set_ylabel('Dwelltime')
    ax.set_title('(B) Dwelltimes')

    # Derived proportions for stacked plots
    total_alive = outcomes_df['total']
    df = pd.DataFrame({
        'years': outcomes_df['years'],
        'prob_clearance': outcomes_df['cleared'] / total_alive * 100,
        'prob_persist': outcomes_df['persisted'] / total_alive * 100,
        'prob_progressed': (outcomes_df['progressed'] + outcomes_df['cancer'] + outcomes_df['dead']) / total_alive * 100,
    })
    total_persisted = outcomes_df['total'] - outcomes_df['cleared']
    df2 = pd.DataFrame({
        'years': outcomes_df['years'],
        'prob_persisted': outcomes_df['persisted'] / total_persisted * 100,
        'prob_progressed': outcomes_df['progressed'] / total_persisted * 100,
        'prob_cancer': (outcomes_df['cancer'] + outcomes_df['dead']) / total_persisted * 100,
    })

    # (C) All outcomes
    ax = axes[2]
    dt_years = df['years'][1] - df['years'][0]
    end_ind = int(1 / dt_years) * 10
    bottom = np.zeros(end_ind)
    for ln, (layer, label) in enumerate([
        ('prob_clearance', 'Cleared'),
        ('prob_persist', 'Infection w/o HSIL'),
        ('prob_progressed', 'HSIL'),
    ]):
        ax.fill_between(df['years'][:end_ind], bottom, bottom + df[layer][:end_ind],
                        color=stacked_colors[ln], label=label)
        bottom += df[layer][:end_ind]
    ax.legend(loc='lower right')
    ax.set_title('(C) All outcomes')
    ax.set_xlabel('Time since infection')

    # (D) Conditional on persistence
    ax = axes[3]
    end_ind = int(1 / dt_years) * 40
    bottom = np.zeros(end_ind)
    for ln, (layer, label) in enumerate([
        ('prob_persisted', 'Infection w/o HSIL'),
        ('prob_progressed', 'HSIL'),
        ('prob_cancer', 'Cancer'),
    ]):
        ax.fill_between(df2['years'][:end_ind], bottom, bottom + df2[layer][:end_ind],
                        color=stacked_colors[ln + 1], label=label)
        bottom += df2[layer][:end_ind]
    ax.legend(loc='lower left')
    ax.set_title('(D) Outcomes conditional on persistence')
    ax.set_xlabel('Time since infection')

    fig.tight_layout()
    sc.savefig(outpath)
    if do_show:
        plt.show()


# %% Run as a script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run the sim and save CSVs (heavy, VM-side)')
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline',
                        help='Dir with plot-ready CSVs (for plot mode only)')
    parser.add_argument('--figpath', default='figures/fig2.png')
    args = parser.parse_args()

    if args.run_sim:
        sim = rs.run_sim(
            calib_pars=None,
            analyzers=[ut.outcomes_by_year(), ut.age_causal(), ut.dwelltime_by_genotype()],
            n_agents=1e6,
            do_save=False,
        )
        save_fig2_data(sim, resfolder='results')
        print('Saved fig2 CSVs to results/ (copy to a versioned baseline dir to commit)')
    else:
        ac_summary = pd.read_csv(f'{args.resfolder}/fig2_age_causal_summary.csv')
        dw_summary = pd.read_csv(f'{args.resfolder}/fig2_dwelltime_summary.csv')
        outcomes_df = pd.read_csv(f'{args.resfolder}/fig2_outcomes.csv')
        plot_fig2(ac_summary, dw_summary, outcomes_df, outpath=args.figpath)
        print('Done.')
