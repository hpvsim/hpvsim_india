'''
Define an HPVsim simulation for Ethiopia
'''

# Standard imports
# Additions to handle numpy multithreading
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import numpy as np
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import behavior_inputs as bi
import utils as ut

# %% Settings and filepaths

# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Run settings  
n_trials    = [10000, 2][debug]  # How many trials to run for calibration
n_workers   = [40, 1][debug]    # How many cores to use
storage     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations

# Save settings
do_save = True
save_plots = True


# %% Simulation creation functions
def make_sim(calib_pars=None, analyzers=[], debug=0, datafile=None, seed=1):
    ''' Define parameters, analyzers, and interventions for the simulation -- not the sim itself '''

    pars = dict(
        n_agents=[50e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        beta=0.28,
        start=[1960, 1980][debug],
        end=2020,
        genotypes=[16, 18, 'hi5', 'ohr'],
        location='india',
        debut=dict(f=dict(dist='lognormal', par1=15., par2=2.),
                   m=dict(dist='lognormal', par1=20., par2=2.)),
        layer_probs=bi.layer_probs,
        m_partners=bi.m_partners,
        f_partners=bi.f_partners,
        f_cross_layer=0.025,
        m_cross_layer=0.25,
        ms_agent_ratio=100,
        verbose=0.0,
    )

    # If calibration parameters have been supplied, use them here
    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    # Create the sim
    sim = hpv.Sim(pars=pars, datafile=datafile, analyzers=analyzers, rand_seed=seed)

    return sim


# %% Simulation running functions
def run_sim(calib_pars=None, analyzers=None, debug=0, datafile=None, seed=1, verbose=.1, do_shrink=True, do_save=False):
    # Make sim
    sim = make_sim(
        debug=debug,
        seed=seed,
        datafile=datafile,
        analyzers=analyzers,
        calib_pars=calib_pars
    )
    sim.label = f'Sim--{seed}'

    # Run
    sim['verbose'] = verbose
    sim.run()
    if do_shrink:
        sim.shrink()

    # Optinally save
    if do_save:
        sim.save(f'results/india.sim')

    return sim


def run_calib(n_trials=None, n_workers=None, do_save=True, filestem=''):

    sim = make_sim()
    datafiles = [
        f'data/india_cancer_cases.csv',
        f'data/india_cin_types.csv',
        f'data/india_cancer_types.csv',
    ]

    # Define the calibration parameters
    genotype_pars = dict(
        hpv16=dict(
            cancer_fn=dict(transform_prob=[2e-3, 1e-3, 3e-3, 2e-4]),
            cin_fn=dict(k=[.3, .2, .4, 0.01]),
            dur_cin=dict(par1=[5, 4, 6, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        hpv18=dict(
            cancer_fn=dict(transform_prob=[2e-3, 1e-3, 3e-3, 2e-4]),
            cin_fn=dict(k=[.25, .15, .35, 0.01]),
            dur_cin=dict(par1=[5, 4, 6, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        hi5=dict(
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
            dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        ohr=dict(
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
            dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
        ),
    )

    calib = hpv.Calibration(sim, calib_pars=None, genotype_pars=genotype_pars,
                            name=f'india_calib',
                            datafiles=datafiles,
                            total_trials=n_trials, n_workers=n_workers,
                            storage=storage
                            )
    calib.calibrate()
    filename = f'india_calib{filestem}'
    if do_save:
        sc.saveobj(f'results/{filename}.obj', calib)

    print(f'Best pars are {calib.best_pars}')

    return sim, calib


def plot_calib(which_pars=0, save_pars=True, filestem=''):
    filename = f'india_calib{filestem}'
    calib = sc.load(f'results/{filename}.obj')

    sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
    sc.options(font='Libertinus Sans')
    fig = calib.plot(res_to_plot=200, plot_type='sns.boxplot', do_save=False)
    fig.tight_layout()
    fig.savefig(f'figures/{filename}.png')

    if save_pars:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
        trial_pars = sc.autolist()
        for i in range(100):
            trial_pars += calib.trial_pars_to_sim_pars(which_pars=i)
        sc.save(f'results//india_pars{filestem}.obj', calib_pars)
        sc.save(f'results/india_pars{filestem}_all.obj', trial_pars)

    return calib


# %% Run as a script
if __name__ == '__main__':

    # List of what to run
    to_run = [
        # 'run_sim',
        'run_calib',
        # 'plot_calib'
    ]

    T = sc.timer()  # Start a timer

    if 'run_sim' in to_run:
        # calib_pars = sc.loadobj('results/india_pars.obj')  # Load parameters from a previous calibration
        sim = run_sim(calib_pars=None, do_shrink=False)  # Run the simulation
        sim.plot()  # Plot the simulation

    if 'run_calib' in to_run:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers, filestem='', do_save=True)

    if 'plot_calib' in to_run:
        calib = plot_calib(save_pars=True, filestem='')

    T.toc('Done')  # Print out how long the run took
