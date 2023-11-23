"""
Plot degree distribution of casual partners
"""
import hpvsim as hpv
import pylab as pl
import utils as ut
import sciris as sc
import numpy as np
import run_sim as rs


# %% Functions
def plot_degree(partners):

    ut.set_font(size=12)
    fig, axes = pl.subplots(1,2, figsize=(9, 5), layout="tight")
    axes = axes.flatten()

    bins = np.concatenate([np.arange(21),[100]]) #np.array([0, 1, 2, 3, 45, 20, 100])

    for ai,sex in enumerate(['f', 'm']):
        counts, bins = np.histogram(partners[sex], bins=bins)
        total = sum(counts)
        counts = counts/total

        axes[ai].bar(bins[:-1], counts)
        axes[ai].set_xlabel(f'Number of lifetime casual partners')
        axes[ai].set_title(f'Distribution of casual partners, {sex}')
        axes[ai].set_ylim([0, 1])
        stats = f"Mean: {np.mean(partners[sex]):.1f}\n"
        stats += f"Median: {np.median(partners[sex]):.1f}\n"
        stats += f"Std: {np.std(partners[sex]):.1f}\n"
        stats += f"%>20: {np.count_nonzero(partners[sex]>=20)/total*100:.2f}\n"
        axes[ai].text(15, 0.5, stats)


    pl.savefig(f"figures/india_degree.png", dpi=300)
    pl.show()

    return

# %% Run as a script
if __name__ == '__main__': 

    do_run = False

    if do_run:
        sim = rs.run_sim(do_shrink=False)
        f_conds = sim.people.is_female * sim.people.alive * sim.people.level0 * sim.people.is_active
        m_conds = sim.people.is_male * sim.people.alive * sim.people.level0 * sim.people.is_active
        partners = {
            'f': sim.people.n_rships[1, f_conds],
            'm': sim.people.n_rships[1, m_conds],
        }
        sc.saveobj('results/partners.obj', partners)
    else:
        partners = sc.loadobj('results/partners.obj')

    plot_degree(partners)

    print('Done.') 
