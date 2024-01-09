"""
Plot implied natural history.
"""
import pylab as pl
import pandas as pd
import numpy as np
import sciris as sc
import utils as ut
import seaborn as sns
import run_sim as rs

do_show = True

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
    ac_colors = sc.gridcolors(3)
    stacked_colors = [
        '#2db5a7',  # cleared
        '#eddc42',  # persisted
        '#e67f2c',  # progressed
        '#871a6c',  # cancer
    ]
    dw_colors = stacked_colors[1:3]

    fig, axes = pl.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()

    ######################################################
    # Top row: dwelltimes and age of causal
    ######################################################
    ac_df = get_age_causal_df(sim)
    ax = axes[0]
    sns.boxplot(
        x="Health event", y="Age", data=ac_df, ax=ax,
        showfliers=False, palette=ac_colors,
    )
    ax.set_title(f'(A) Age distribution of key health events')
    ax.set_xlabel('')
    ax.set_ylim([0, 100])

    # Panel B
    dw_df = get_dwelltime_df(sim=sim)
    ax = axes[1]
    sns.boxplot(data=dw_df, x="state", y="dwelltime",  ax=ax,
                showfliers=False, palette=dw_colors)
    ax.legend([], [], frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Dwelltime")
    ax.set_title('(B) Dwelltimes')


    ######################################################
    # Bottom row: stacked plots
    ######################################################
    a = sim.get_analyzer('outcomes_by_year')
    res = a.results
    years = a.durations

    df = pd.DataFrame()
    total_alive = res["total"]
    df["years"] = years
    df["prob_clearance"] = (res["cleared"]) / total_alive * 100
    df["prob_persist"] = (res["persisted"]) / total_alive * 100
    df["prob_progressed"] = (res["progressed"] + res["cancer"] + res["dead"]) / total_alive * 100

    df2 = pd.DataFrame()
    total_persisted = res["total"] - res["cleared"]
    df2["years"] = years
    df2["prob_persisted"] = (res["persisted"]) / total_persisted * 100
    df2["prob_progressed"] = (res["progressed"]) / total_persisted * 100
    df2["prob_cancer"] = (res["cancer"] + res["dead"]) / total_persisted * 100

    ####################
    # Make figure, set fonts and colors
    ####################

    # Panel C, all outcomes
    ax = axes[2]
    n_years = 10
    end_ind = int(1/(a.durations[1]-a.durations[0]))*n_years
    bottom = np.zeros(len(df["years"][:end_ind]))
    layers = [
        "prob_clearance",
        "prob_persist",
        "prob_progressed",
    ]
    labels = ["Cleared", "Infection w/o HSIL", "HSIL"]
    for ln, layer in enumerate(layers):
        ax.fill_between(
            df["years"][:end_ind], bottom, bottom + df[layer][:end_ind], color=stacked_colors[ln], label=labels[ln]
        )
        bottom += df[layer][:end_ind]
    ax.legend(loc="lower right")
    ax.set_title('(C) All outcomes')
    ax.set_xlabel("Time since infection")

    # Panel D, conditional on being alive and not cleared
    ax = axes[3]
    n_years = 40
    end_ind = int(1/(a.durations[1]-a.durations[0]))*n_years
    bottom = np.zeros(len(df2["years"]))
    layers = ["prob_persisted", "prob_progressed", "prob_cancer"]
    labels = ["Infection w/o HSIL", "HSIL", "Cancer"]
    for ln, layer in enumerate(layers):
        ax.fill_between(
            df2["years"][:end_ind],
            bottom[:end_ind],
            bottom[:end_ind] + df2[layer][:end_ind],
            color=stacked_colors[ln+1],
            label=labels[ln],
        )
        bottom += df2[layer]
    ax.legend(loc="lower left")
    ax.set_title('(D) Outcomes conditional on persistence')
    ax.set_xlabel("Time since infection")

    fig.tight_layout()
    sc.savefig(f"figures/fig2.png")

    if do_show:
        pl.show()

    return


# %% Run as a script
if __name__ == '__main__':

    location = 'india'
    make_sim = True # Needed since has a different analyzer
    if make_sim:
        # calib_pars = sc.loadobj('results/india_pars.obj')  # Load parameters from a previous calibration
        sim = rs.run_sim(
            calib_pars=None,
            analyzers=[ut.outcomes_by_year(), ut.age_causal(), ut.dwelltime_by_genotype()],
            n_agents=1e6,
            do_save=False)  # Run the simulation
    else:
        sim = sc.loadobj(f'results/{location}.sim')

    plot_fig2(sim)

    print('Done.')
