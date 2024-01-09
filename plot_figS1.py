"""
Plot sexual behavior data
"""

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns
from scipy.stats import norm, lognorm
import hpvsim as hpv
import hpvsim.utils as hpu

# Imports from this repository
import utils as ut
import run_sim as rs

do_show = True


def plot_sb(dist_type='lognormal'):
    '''
    Create plots of sexual behavior inputs and outputs
    '''

    ut.set_font(15)
    fig = pl.figure(layout="tight", figsize=(11, 10))

    # Make 2 rows, with 3 panels in the top row and 2 in the bottom
    gs0 = fig.add_gridspec(2, 1)
    gs00 = gs0[0].subgridspec(1, 3)
    gs01 = gs0[1].subgridspec(1, 2)
    ms = 80

    # Panel A: debut
    data_countries, dff, df2, rvs = ut.read_debut_data(dist_type=dist_type)
    alldf = sc.loadobj(f'results/model_sb_AFS.obj')
    sex = 'Women'
    sk = 'f'
    dfw = dff[sex]

    ax = fig.add_subplot(gs00[0])
    dfplot = dfw.loc[(dfw["AgeStr"]!=f'{sex} never')&(dfw["AgeStr"]!=f'{sex} 60')]
    dfmed = df2

    rv = rvs[sex]['India']
    xx = np.arange(12,30.1,0.1)
    xxx = np.arange(12,31,1)

    for cohort in alldf["cohort"].unique():
        modely = np.array(alldf.loc[(alldf["cohort"]==cohort)][f'model_prop_{sk}'])
        ax.plot(xxx, modely*100, 'b-', lw=1, alpha=0.3)

    sns.scatterplot(ax=ax, data=dfplot, x="Age", y="Percentage", marker='d', s=ms, color='k')
    sns.scatterplot(ax=ax, data=dfmed, x=f"{sex} median", y="y", marker='d', s=ms, color='k')
    ax.plot(xx, rv.cdf(xx)*100, 'k--', lw=2)

    ax.set_ylabel('Share')
    ax.set_xlabel('Age')
    ax.set_title('(A) Share of females who\n are sexually active')

    # Panel B: proportion married
    dfraw = pd.read_csv('data/prop_married.csv')
    df = dfraw.melt(id_vars=['Country', 'Survey'], value_name='Percentage', var_name='AgeRange')
    modeldf = sc.loadobj(f'results/model_sb_prop_married.obj')
    modeldf.reset_index()

    colors = sc.gridcolors(1)

    ax = fig.add_subplot(gs00[1])

    sns.scatterplot(ax=ax, data=df, x="AgeRange", y="Percentage")

    modeldf['val'] = modeldf['val'].apply(lambda x: x * 100)
    sns.boxplot(data=modeldf, x="age", y="val", color=colors[0], ax=ax)
    ax.set_ylabel('Share')
    ax.set_xlabel('Age')
    ax.set_title('(B) Share of females\n who are married')

    # Panel C of age diffs
    agediffs = sc.loadobj(f'results/model_age_diffs.obj')
    ax = fig.add_subplot(gs00[2])

    # Plot model
    dfplot_m = agediffs
    sns.kdeplot(data=dfplot_m, color=colors[0], ax=ax)
    ax.legend([], [], frameon=False)
    ax.set_xlim([-10, 30])
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_ylabel('Share')
    ax.set_xlabel('Male age - female age')
    ax.set_title('(C) Age differences\n between partners')

    ########################################
    # Bottom row - degree distribution
    ########################################
    bins = np.concatenate([np.arange(21),[100]]) #np.array([0, 1, 2, 3, 45, 20, 100])
    partners = sc.loadobj('results/partners.obj')
    axlabels = ['D', 'E']

    for ai,slabel in enumerate(['females', 'males']):
        sex = slabel[0]
        counts, bins = np.histogram(partners[sex], bins=bins)
        total = sum(counts)
        counts = counts/total

        ax = fig.add_subplot(gs01[ai])
        ax.bar(bins[:-1], counts)
        ax.set_xlabel(f'Number of lifetime casual partners')
        ax.set_title(f'({axlabels[ai]}) Distribution of casual partners, {slabel}')
        ax.set_ylim([0, 1])
        stats = f"Mean: {np.mean(partners[sex]):.1f}\n"
        stats += f"Median: {np.median(partners[sex]):.1f}\n"
        stats += f"Std: {np.std(partners[sex]):.1f}\n"
        stats += f"%>20: {np.count_nonzero(partners[sex]>=20)/total*100:.2f}\n"
        ax.text(15, 0.5, stats)


    fig.tight_layout()
    sc.savefig(f"figures/india_behavior.png", dpi=100)
    if do_show:
        pl.show()

    return


#%% Run as a script
if __name__ == '__main__':

    plot_sb()

    print('Done.')
