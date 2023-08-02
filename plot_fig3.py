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


# %% Functions

def plot_nh(sim=None):
    # Make sims
    genotypes = ['hpv16', 'hpv18', 'hi5', 'ohr']
    glabels = ['HPV16', 'HPV18', 'HI5', 'OHR']
    clinical_cutoffs = dict(cin1=0.33, cin2=0.676, cin3=0.8)

    dur_episomal = sc.autolist()
    transform_probs = sc.autolist()
    sev_fns = sc.autolist()
    dur_precin = sc.autolist()
    for gi, genotype in enumerate(genotypes):
        dur_precin += sim['genotype_pars'][genotype]['dur_precin']
        dur_episomal += sim['genotype_pars'][genotype]['dur_episomal']
        transform_probs += sim['genotype_pars'][genotype]['transform_prob']
        sev_fns += sim['genotype_pars'][genotype]['sev_fn']
        # sims[location] = sim
    sev_dist = sim['sev_dist']

    ####################
    # Make figure, set fonts and colors
    ####################
    ut.set_font(size=16)
    colors = sc.gridcolors(len(genotypes))
    fig, axes = pl.subplots(2, 3, figsize=(11, 9))
    axes = axes.flatten()
    cmap = pl.cm.Oranges([0.25, 0.5, 0.75, 1])

    ####################
    # Make plots
    ####################
    dt = 0.25
    max_x = 30
    thisx = np.arange(dt, max_x+dt, dt)
    n_samples = 10
    # Durations and severity of dysplasia
    for gi, gtype in enumerate(genotypes):

        # Panel A: durations
        sigma, scale = ut.lognorm_params(dur_episomal[gi]['par1'], dur_episomal[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        axes[0].plot(thisx, rv.pdf(thisx), color=colors[gi], lw=2, label=glabels[gi])

        # Panel C: dysplasia
        dysp = np.zeros_like(thisx)
        dysp_start_ind = sc.findnearest(thisx, dur_precin[gi]['par1'])
        if dysp_start_ind > 0:
            dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], pars=sev_fns[gi])
        else:
            dysp[:] = hppar.compute_severity(thisx[:], pars=sev_fns[gi])
        axes[1].plot(thisx, dysp, color=colors[gi], lw=2, label=gtype.upper())

        # Panel B: dysplasia
        dysp = np.zeros_like(thisx)
        dysp_start_ind = sc.findnearest(thisx, dur_precin[gi]['par1'])
        if dysp_start_ind > 0:
            dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], pars=sev_fns[gi])
        else:
            dysp[:] = hppar.compute_severity(thisx[:], pars=sev_fns[gi])
        axes[1].plot(thisx, dysp, color=colors[gi], lw=2, label=gtype.upper())

        # Panel C: transform prob
        cum_dysp = np.zeros_like(thisx)
        dysp_int = hppar.compute_severity_integral(thisx, pars=sev_fns[gi])

        if dysp_start_ind > 0:
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            cum_dysp[:] = dysp_int[:]

        tp_array = hpu.transform_prob(transform_probs[gi], cum_dysp)
        axes[2].plot(thisx, tp_array, color=colors[gi], lw=2, label=gtype.upper())

        # # Add rel_sev samples to B and C plots
        # for smpl in range(n_samples):
        #     rel_sev = hpu.sample(**sev_dist)
        #     dysp_start = hpu.sample(**dur_precin[gi])
        #     rel_dysp = np.zeros_like(thisx)
        #     rel_cum_dysp = np.zeros_like(thisx)
        #     dysp_start_ind = sc.findnearest(thisx, dysp_start)
        #     rel_dysp_int = hppar.compute_severity_integral(thisx, rel_sev=rel_sev, pars=sev_fns[gi])
        #
        #     if dysp_start_ind>0:
        #         rel_dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], rel_sev=rel_sev, pars=sev_fns[gi])
        #         rel_cum_dysp[dysp_start_ind:] = rel_dysp_int[:-dysp_start_ind]
        #     elif dysp_start_ind==0:
        #         rel_dysp[:] = hppar.compute_severity(thisx[:], rel_sev=rel_sev, pars=sev_fns[gi])
        #         rel_cum_dysp[:] = rel_dysp_int[:]
        #     axes[1].plot(thisx, rel_dysp, color=colors[gi], lw=1, alpha=0.5)
        #
        #     rel_tp_array = hpu.transform_prob(transform_probs[gi], rel_cum_dysp)
        #     axes[2].plot(thisx, rel_tp_array, color=colors[gi], lw=1, alpha=0.5)

    axes[0].set_ylabel("")
    axes[0].grid()
    axes[0].set_xlabel("Duration of infection (years)")
    axes[0].set_title("Distribution of\n infection duration")
    axes[0].legend(frameon=False)

    axes[1].set_ylabel("Severity of infection")
    axes[1].set_xlabel("Duration of infection (years)")
    axes[1].set_title("Duration to severity\nfunction")
    axes[1].set_ylim([0,1])
    axes[1].grid()

    axes[1].axhline(y=clinical_cutoffs['cin1'], ls=':', c='k')
    axes[1].axhline(y=clinical_cutoffs['cin2'], ls=':', c='k')
    axes[1].axhspan(0, clinical_cutoffs['cin1'], color=cmap[0], alpha=.4)
    axes[1].axhspan(clinical_cutoffs['cin1'], clinical_cutoffs['cin2'], color=cmap[1], alpha=.4)
    axes[1].axhspan(clinical_cutoffs['cin2'], 1.0, color=cmap[2], alpha=.4)
    axes[1].text(-0.3, 0.15, 'CIN1', rotation=90)
    axes[1].text(-0.3, 0.48, 'CIN2', rotation=90)
    axes[1].text(-0.3, 0.8, 'CIN3', rotation=90)

    axes[2].grid()
    # axes[2].set_ylabel("Probability of transformation")
    axes[2].set_xlabel("Duration of infection (years)")
    axes[2].set_title("Probability of cancer\n within X years")


    dt = 0.25
    max_x = 30
    x = np.arange(dt, max_x + dt, dt)
    annual_x = np.arange(1, 11, 1)
    width = 0.2  # the width of the bars
    multiplier = 0

    # Panel A: clearance rates
    for gi, genotype in enumerate(genotypes):
        offset = width * multiplier
        sigma, scale = ut.lognorm_params(dur_episomal[gi]['par1'], dur_episomal[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        axes[3].bar(annual_x + offset - width / 2, rv.cdf(annual_x), color=colors[gi], lw=2, label=genotype.upper(),
                    width=width)
        multiplier += 1
    axes[3].set_title("Proportion clearing\n within X years")
    axes[3].set_xticks(annual_x)
    axes[3].set_ylabel("Probability")
    axes[3].set_xlabel("Years")

    # # Panel B: transform prob
    # for gi, genotype in enumerate(genotypes):
    #     cum_dysp = np.zeros_like(x)
    #     dysp_int = hppar.compute_severity_integral(x, pars=sev_fns[gi])
    #     dysp_start_ind = sc.findnearest(x, dur_precins[gi]['par1'])
    #     if dysp_start_ind > 0:
    #         cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
    #     else:
    #         cum_dysp[:] = dysp_int[:]
    #     tp_array = hpu.transform_prob(transform_probs[gi], cum_dysp)
    #     axes[1].plot(x, tp_array, color=colors[gi], lw=2, label=genotype.upper())
    # axes[1].set_title("Probability of cancer\n within X years")
    # axes[1].set_ylabel("Probability")
    # axes[1].set_xlabel("Years")
    # axes[1].legend()

    # Panel B: CIN dwelltime
    a = sim.get_analyzer('dwelltime_by_genotype')
    dd = {}
    dd['dwelltime'] = sc.autolist()
    dd['genotype'] = sc.autolist()
    dd['state'] = sc.autolist()
    for cin in ['cin1', 'cin2', 'cin3']:
        dt = a.dwelltime[cin]
        data = dt[0]+dt[1]+dt[2]+dt[3]
        labels = ['HPV16']*len(dt[0]) + ['HPV18']*len(dt[1]) + ['HI5']*len(dt[2]) + ['OHR']*len(dt[3])
        dd['dwelltime'] += data
        dd['genotype'] += labels
        dd['state'] += [cin.upper()]*len(labels)
    df = pd.DataFrame(dd)
    # sns.violinplot(data=df, x="state", y="dwelltime", hue="genotype", ax=axes[1], cut=0)
    sns.boxplot(data=df, x="state", y="dwelltime", hue="genotype", ax=axes[4], showfliers=False, palette=colors)
    # handles, labels = axes[1].get_legend_handles_labels()
    # axes[4].legend(handles, labels, frameon=False)
    axes[4].legend([], [], frameon=False)
    axes[4].set_xlabel("")
    axes[4].set_ylabel("Dwelltime")
    axes[4].set_title('Dwelltimes from\n infection to CIN grades')


    # sns.violinplot(data=dd, x="genotype", y="dwelltime", ax=axes[2], palette=colors)
    # axes[2].set_xlabel('')
    # axes[2].set_ylabel('')
    # axes[2].set_ylabel("Years")
    # axes[2].set_title('Total dwelltime\n from infection to cancer')


    # Panel C: total dwelltime
    dd = pd.DataFrame()
    dw = sc.autolist()
    gen = sc.autolist()
    for gi, genotype in enumerate(genotypes):
        a = sim.get_analyzer('dwelltime_by_genotype')
        dw += a.dwelltime['total'][gi]
        gen += [genotype.upper()] * len(a.dwelltime['total'][gi])
    dd['genotype'] = gen
    dd['dwelltime'] = dw
    sns.violinplot(data=dd, x="genotype", y="dwelltime", ax=axes[5], palette=colors)
    axes[5].set_xlabel('')
    axes[5].set_ylabel('')
    axes[5].set_ylabel("Years")
    axes[5].set_title('Total dwelltime\n from infection to cancer')

    fig.tight_layout()
    fs=24
    pl.figtext(0.04, 0.955, 'A', fontweight='bold', fontsize=fs)
    pl.figtext(0.36, 0.955, 'B', fontweight='bold', fontsize=fs)
    pl.figtext(0.7, 0.955, 'C', fontweight='bold', fontsize=fs)
    pl.figtext(0.04, 0.47, 'D', fontweight='bold', fontsize=fs)
    pl.figtext(0.36, 0.47, 'E', fontweight='bold', fontsize=fs)
    pl.figtext(0.7, 0.47, 'F', fontweight='bold', fontsize=fs)

    pl.savefig(f"figures/fig3.png", dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    location = 'india'
    make_sim = True
    if make_sim:
        calib_pars = sc.loadobj('results/india_pars.obj')  # Load parameters from a previous calibration
        sim = rs.run_sim(calib_pars=None, analyzers=ut.dwelltime_by_genotype(), do_save=True)  # Run the simulation
    else:
        sim = sc.loadobj(f'results/{location}.sim')

    plot_nh(sim)

    print('Done.')
