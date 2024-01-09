"""
Plot India metrics from HPVsim and Globocan
"""
import hpvsim as hpv
import pylab as pl
import pandas as pd
import numpy as np
import sciris as sc
import utils as ut
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter


do_show = True

# %% Functions
def plot_india_burden():

    ut.set_font(size=16)
    res = sc.loadobj(f'results/india_msim.obj')

    colors = sc.gridcolors(2)

    globocan = pd.read_csv(f'data/globocan_asr_cancer_incidence.csv')

    fig, ax = pl.subplots(1, 1, figsize=(10, 5))
    start_year = 2005
    ind = sc.findinds(res['year'], start_year)[0]
    years = res['year'][ind:]

    ax.plot(years, res['asr_cancer_incidence'].values[ind:], color=colors[0], label=f'HPVsim', marker='o')
    ax.fill_between(years, res['asr_cancer_incidence'].low[ind:],
                            res['asr_cancer_incidence'].high[ind:], color=colors[0], alpha=0.3)

    gc_ind = sc.findinds(globocan['year'], start_year)[0]
    ax.plot(globocan['year'][gc_ind:], globocan['asr_cancer_incidence'][gc_ind:], marker='s', color=colors[1], label='Globocan')
    ax.set_ylabel('Age-standardized cancer incidence (per 100k)')
    ax.legend()
    ax.set_ylim([0,25])
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    # sc.SIticks(ax)
    fig.tight_layout()
    fig_name = f'figures/figS3.png'
    sc.savefig(fig_name, dpi=100)
    if do_show:
        pl.show()

    return

# %% Run as a script
if __name__ == '__main__':


    plot_india_burden()

    print('Done.') 
