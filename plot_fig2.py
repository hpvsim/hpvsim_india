"""
Plot implied natural history.
"""
import hpvsim as hpv
import hpvsim.utils as hpu
import hpvsim.parameters as hppar
import pylab as pl
import pandas as pd
from scipy.stats import lognorm, norm
import numpy as np
import sciris as sc
import utils as ut
import seaborn as sns

import run_sim as rs
import math

# %% Functions
def get_age_causal_df(sim=None):
    """
    Make age causal dataframe
    """
    dt_res = sim.get_analyzer('age_causal')
    dt_dfs = sc.autolist()

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_causal)[np.array(dt_res.age_causal)<50]
    dt_df['Health event'] = 'Causal infection'
    dt_dfs += dt_df

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_cin)[np.array(dt_res.age_causal)<65]
    dt_df['Health event'] = 'CIN2+'
    dt_dfs += dt_df

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_cancer)[np.array(dt_res.age_causal)<90]
    dt_df['Health event'] = 'Cancer'
    dt_dfs += dt_df

    age_causal_df = pd.concat(dt_dfs)

    return age_causal_df


def get_dwelltime_df(sim):
    a = sim.get_analyzer('dwelltime_by_genotype')
    dd = {}
    dd['dwelltime'] = sc.autolist()
    dd['genotype'] = sc.autolist()
    dd['state'] = sc.autolist()
    for cin in ['precin', 'cin']:
        dt = a.dwelltime[cin]
        data = dt[0]+dt[1]+dt[2]+dt[3]
        labels = ['HPV16']*len(dt[0]) + ['HPV18']*len(dt[1]) + ['HI5']*len(dt[2])+ ['OHR']*len(dt[3])
        dd['dwelltime'] += data
        dd['genotype'] += labels
        dd['state'] += [cin.upper()]*len(labels)
    df = pd.DataFrame(dd)
    return df


def plot_fig2(sim):

    ut.set_font(size=16)
    ac_color = '#a8327b'
    dw_color = '#027d23'
    colors = sc.gridcolors(5)
    fig, axes = pl.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()

    ######################################################
    # Top row: dwelltimes and age of causal
    ######################################################
    ac_df = get_age_causal_df(sim)
    ax = axes[0]
    sns.boxplot(
        x="Health event", y="Age", data=ac_df, ax=ax,
        showfliers=False, color=ac_color,
    )
    ax.set_title(f'(A) Age distribution of key health events')
    ax.set_xlabel('')
    ax.set_ylim([0, 100])

    # Panel B
    dw_df = get_dwelltime_df(sim=sim)
    ax = axes[1]
    sns.boxplot(data=dw_df, x="state", y="dwelltime",  ax=ax,
                showfliers=False, color=dw_color)
    ax.legend([], [], frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Dwelltime")
    ax.set_title('(B) Dwelltimes')


    ######################################################
    # Bottom row: stacked plots
    ######################################################

    # # Jamie's version
    # cum_dist = sim.get_analyzer('cum_dist')
    # durs_to_cancer, counts_to_cancer = np.unique([math.ceil(elem) for elem in cum_dist.dur_to_cancer], return_counts=True)
    # durs_to_cin, counts_to_cin = np.unique([math.ceil(elem) for elem in cum_dist.dur_to_cin], return_counts=True)
    # durs_to_clearance, counts_to_clearance = np.unique([math.ceil(elem) for elem in cum_dist.dur_to_clearance], return_counts=True)
    #
    # df = pd.DataFrame()
    # df['years'] = np.arange(0, 30)
    # durs = np.zeros(30)
    # durs_subset = durs_to_clearance[durs_to_clearance < 30]
    # durs[[int(elem) for elem in durs_subset]] = counts_to_clearance[:len(durs_subset)]
    # df['n_cleared'] = durs
    # df['prob_clearance'] = 100 * np.cumsum(df['n_cleared']) / cum_dist.total_infections
    #
    # durs_subset = durs_to_cin[durs_to_cin < 30]
    # durs[[int(elem) for elem in durs_subset]] = counts_to_cin[:len(durs_subset)]
    # df['n_cin'] = durs
    # df['prob_cin'] = 100 * np.cumsum(df['n_cin']) / cum_dist.total_infections
    #
    # durs_subset = durs_to_cancer[durs_to_cancer < 30]
    # durs[[int(elem) for elem in durs_subset]] = counts_to_cancer[:len(durs_subset)]
    # df['n_cancer'] = durs
    # df['prob_cancer'] = 100 * np.cumsum(df['n_cancer']) / cum_dist.total_infections
    #
    # ax = axes[2]
    # ax.fill_between(df['years'], np.zeros(len(df['years'])), df['prob_clearance'], color=colors[0], label='Cleared')
    # ax.fill_between(df['years'], df['prob_clearance'], 100 - df['prob_cin'], color=colors[1], label='Persisted')
    # ax.fill_between(df['years'], 100 - df['prob_cin'], 100 - df['prob_cancer'], color=colors[2], label='CIN2+')
    # ax.fill_between(df['years'], 100 - df['prob_cancer'], 100 * np.ones(len(df['years'])), color=colors[3], label='Cancer')
    # ax.legend(loc='lower right')
    # ax.set_xlabel('Time since infection')


    a = sim.get_analyzer('outcomes_by_year')
    res = a.results
    years = a.durations

    df = pd.DataFrame()
    total_alive = res["total"] - res["dead"]
    df["years"] = years
    df["prob_clearance"] = (res["cleared"]) / total_alive * 100
    df["prob_persist"] = (res["persisted"]) / total_alive * 100
    df["prob_progressed"] = (res["progressed"] + res["cancer"] + res["dead"]) / total_alive * 100

    df2 = pd.DataFrame()
    total_persisted = res["total"] - res["cleared"]
    df2["years"] = years
    df2["prob_persisted"] = (res["persisted"] + res["progressed"]) / total_persisted * 100
    df2["prob_cancer"] = (res["cancer"] + res["dead"]) / total_persisted * 100

    ####################
    # Make figure, set fonts and colors
    ####################

    # Panel C, all outcomes
    ax = axes[2]
    bottom = np.zeros(len(df["years"][0:10]))
    layers = [
        "prob_clearance",
        "prob_persist",
        "prob_progressed",
    ]
    labels = ["Cleared", "Persistent Infection", "CIN2+"]
    for ln, layer in enumerate(layers):
        ax.fill_between(
            df["years"][0:10], bottom, bottom + df[layer][0:10], color=colors[ln], label=labels[ln]
        )
        bottom += df[layer][0:10]
    ax.legend(loc="lower right")
    ax.set_title('(C) Short-term outcomes')
    ax.set_xlabel("Time since infection")

    # Panel D, conditional on being alive and not cleared
    ax = axes[3]
    bottom = np.zeros(len(df2["years"]))
    layers = ["prob_persisted", "prob_cancer"]
    labels = ["CIN2+ regression/persistence", "Cancer"]
    for ln, layer in enumerate(layers):
        ax.fill_between(
            df2["years"],
            bottom,
            bottom + df2[layer],
            color=colors[ln+2],
            label=labels[ln],
        )
        bottom += df2[layer]
    ax.legend(loc="lower right")
    ax.set_title('(D) Long-term outcomes')
    ax.set_xlabel("Time since infection")

    fig.tight_layout()
    fig.savefig(f"figures/fig2.png")
    fig.show()


    return


# %% Run as a script
if __name__ == '__main__':

    location = 'india'
    make_sim = False
    if make_sim:
        # calib_pars = sc.loadobj('results/india_pars.obj')  # Load parameters from a previous calibration
        sim = rs.run_sim(calib_pars=None, analyzers=[ut.outcomes_by_year(), ut.cum_dist(), ut.age_causal(), ut.dwelltime_by_genotype()], do_save=True)  # Run the simulation
    else:
        sim = sc.loadobj(f'results/{location}.sim')

    plot_fig2(sim)

    print('Done.')
