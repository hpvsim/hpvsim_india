"""
Plot sexual behavior data.

Two modes:
  python plot_figS1.py --run-sim   # run sim + extract sexual-behavior artifacts (VM)
  python plot_figS1.py             # plot from saved CSVs (local)
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


def plot_sb(dist_type='lognormal', resfolder='results', outpath='figures/india_behavior.png'):
    ut.set_font(15)
    fig = plt.figure(layout='tight', figsize=(11, 10))
    gs0 = fig.add_gridspec(2, 1)
    gs00 = gs0[0].subgridspec(1, 3)
    gs01 = gs0[1].subgridspec(1, 2)
    ms = 80

    # Panel A: age at first sex
    data_countries, dff, df2, rvs = ut.read_debut_data(dist_type=dist_type)
    alldf = pd.read_csv(f'{resfolder}/model_sb_AFS.csv')
    sex = 'Women'
    sk = 'f'
    dfw = dff[sex]

    ax = fig.add_subplot(gs00[0])
    dfplot = dfw.loc[(dfw['AgeStr'] != f'{sex} never') & (dfw['AgeStr'] != f'{sex} 60')]
    dfmed = df2

    rv = rvs[sex]['India']
    xx = np.arange(12, 30.1, 0.1)
    xxx = np.arange(12, 31, 1)

    for cohort in alldf['cohort'].unique():
        modely = alldf.loc[alldf['cohort'] == cohort, f'model_prop_{sk}'].values
        ax.plot(xxx, modely * 100, 'b-', lw=1, alpha=0.3)

    sns.scatterplot(ax=ax, data=dfplot, x='Age', y='Percentage', marker='d', s=ms, color='k')
    sns.scatterplot(ax=ax, data=dfmed, x=f'{sex} median', y='y', marker='d', s=ms, color='k')
    ax.plot(xx, rv.cdf(xx) * 100, 'k--', lw=2)
    ax.set_ylabel('Share')
    ax.set_xlabel('Age')
    ax.set_title('(A) Share of females who\n are sexually active')

    # Panel B: proportion married
    modeldf = pd.read_csv(f'{resfolder}/model_sb_prop_married.csv')
    modeldf['val'] = modeldf['val'] * 100
    modeldf['AgeRange'] = modeldf['age'].astype(str) + '-' + (modeldf['age'] + 4).astype(str)

    colors = sc.gridcolors(1)
    ax = fig.add_subplot(gs00[1])
    sns.boxplot(data=modeldf, x='AgeRange', y='val', color=colors[0], ax=ax,
                order=sorted(modeldf['AgeRange'].unique(), key=lambda s: int(s.split('-')[0])))
    try:
        dfraw = pd.read_csv('data/prop_married.csv')
        df = dfraw.melt(id_vars=['Country', 'Survey'], value_name='Percentage', var_name='AgeRange')
        df_loc = df[df['Country'] == 'India']
        if len(df_loc):
            sns.scatterplot(ax=ax, data=df_loc, x='AgeRange', y='Percentage', color='k', marker='d', s=ms)
    except FileNotFoundError:
        pass
    ax.set_ylabel('Share')
    ax.set_xlabel('Age')
    ax.set_title('(B) Share of females\n who are married')

    # Panel C: age differences between partners (precomputed kde grid)
    kde_df = pd.read_csv(f'{resfolder}/age_diffs_kde.csv')
    ax = fig.add_subplot(gs00[2])
    ax.plot(kde_df['x'], kde_df['density'], color=colors[0])
    ax.set_xlim([-10, 30])
    ax.set_ylabel('Share')
    ax.set_xlabel('Male age - female age')
    ax.set_title('(C) Age differences\n between partners')

    # Panels D & E: degree distribution (precomputed histogram)
    partners_hist = pd.read_csv(f'{resfolder}/partners_hist.csv')
    axlabels = ['D', 'E']
    for ai, slabel in enumerate(['females', 'males']):
        s = slabel[0]
        sub = partners_hist[partners_hist['sex'] == s].sort_values('bin')
        ax = fig.add_subplot(gs01[ai])
        ax.bar(sub['bin'].values, sub['probability'].values)
        ax.set_xlabel('Number of lifetime casual partners')
        ax.set_title(f'({axlabels[ai]}) Distribution of casual partners, {slabel}')
        ax.set_ylim([0, 1])
        row = sub.iloc[0]
        stats = (
            f'Mean: {row["mean"]:.1f}\n'
            f'Median: {row["median"]:.1f}\n'
            f'Std: {row["std"]:.1f}\n'
            f'%>20: {row["pct_gt_20"]:.2f}\n'
        )
        ax.text(15, 0.5, stats)

    fig.tight_layout()
    sc.savefig(outpath, dpi=100)
    if do_show:
        plt.show()


# %% Run as a script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run the sim and save sexual-behavior CSVs (VM-side)')
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline',
                        help='Dir with plot-ready CSVs (for plot mode only)')
    parser.add_argument('--figpath', default='figures/india_behavior.png')
    args = parser.parse_args()

    if args.run_sim:
        rs.get_sb_from_sims()
        print(f'Saved sexual-behavior CSVs to results/')
    else:
        plot_sb(resfolder=args.resfolder, outpath=args.figpath)
        print('Done.')
