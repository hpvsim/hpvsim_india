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
import pandas as pd
import seaborn as sns
import pylab as pl

# Imports from this repository
import behavior_inputs as bi

# %% Settings and filepaths

# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Run settings
n_trials    = [8000, 2][debug]  # How many trials to run for calibration
n_workers   = [40, 4][debug]    # How many cores to use
storage     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug] # Storage for calibrations

# Save settings
do_save = True
save_plots = True


# %% Simulation creation functions
def make_sim(calib_pars=None, debug=0, datafile=None, seed=1):
    ''' Define parameters, analyzers, and interventions for the simulation -- not the sim itself '''

    pars = dict(
        n_agents=[10e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        end=2020,
        network='default',
        genotypes=[16, 18, 'hi5', 'ohr'],
        location='india',
        debut=dict(f=dict(dist='lognormal', par1=14.8, par2=2.),
                   m=dict(dist='lognormal', par1=17.0, par2=2.)),
        mixing=bi.default_mixing,
        layer_probs=bi.default_layer_probs,
        partners=bi.default_partners,
        init_hpv_dist=dict(hpv16=0.4, hpv18=0.15, hi5=0.15, ohr=0.3),
        init_hpv_prev={
            'age_brackets': np.array([12, 17, 24, 34, 44, 64, 80, 150]),
            'm': np.array([0.0, 0.25, 0.6, 0.25, 0.05, 0.01, 0.0005, 0]),
            'f': np.array([0.0, 0.35, 0.7, 0.25, 0.05, 0.01, 0.0005, 0]),
        },
        ms_agent_ratio=100,
        verbose=0.0,
    )

    # If calibration parameters have been supplied, use them here
    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    # Create the sim
    sim = hpv.Sim(pars=pars, datafile=datafile, rand_seed=seed)

    return sim


# %% Simulation running functions
def run_sim(calib_pars=None, debug=0, datafile=None, seed=1, verbose=.1, do_save=False):
    # Make sim
    sim = make_sim(
        debug=debug,
        seed=seed,
        datafile=datafile,
        calib_pars=calib_pars
    )
    sim.label = f'Sim--{seed}'

    # Run
    sim['verbose'] = verbose
    sim.run()
    sim.shrink()

    # Optinally save
    if do_save:
        sim.save(f'results/india.sim')

    return sim


def run_calib(n_trials=None, n_workers=None, do_save=True, filestem=''):

    sim = make_sim()
    datafiles = [
        f'data/india_hpv_prevalence.csv',
        f'data/india_cancer_cases.csv',
        f'data/india_cin1_types.csv',
        f'data/india_cin3_types.csv',
        f'data/india_cancer_types.csv',
    ]

    # Define the calibration parameters
    calib_pars = dict(
        beta=[0.05, 0.02, 0.5, 0.01],
    )
    genotype_pars = dict(
        hpv16=dict(
            transform_prob=[10e-10, 4e-10, 20e-10, 1e-10],
            sev_fn=dict(k=[0.25, 0.15, 0.4, 0.05]),
        ),
        hpv18=dict(
            transform_prob=[6e-10, 4e-10, 10e-10, 1e-10],
            sev_fn=dict(k=[0.2, 0.1, 0.35, 0.05]),
        ),
        hi5=dict(
                transform_prob=[3e-10, 2e-10, 5e-10, 1e-10],
                sev_fn=dict(k=[0.05, 0.04, 0.2, 0.01]),
            ),
        ohr=dict(
            transform_prob=[3e-10, 2e-10, 5e-10, 1e-10],
            sev_fn=dict(k=[0.05, 0.04, 0.2, 0.01]),
        ),
    )

    calib = hpv.Calibration(sim, calib_pars=calib_pars, genotype_pars=genotype_pars,
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
        'run_sim',
        # 'run_calib',
        # 'plot_calib'
    ]

    T = sc.timer()  # Start a timer

    if 'run_sim' in to_run:
        # calib_pars = sc.loadobj('results/india_pars.obj')  # Load parameters from a previous calibration
        sim = run_sim(calib_pars=None)  # Run the simulation
        sim.plot()  # Plot the simulation

    if 'run_calib' in to_run:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers, do_save=do_save, filestem='')

    if 'plot_calib' in to_run:
        calib = plot_calib(save_pars=True, filestem='')

    T.toc('Done')  # Print out how long the run took
