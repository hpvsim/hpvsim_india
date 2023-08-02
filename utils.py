"""
Utilities
"""

# Imports
import sciris as sc
import numpy as np
import hpvsim as hpv


def set_font(size=None, font='Libertinus Sans'):
    """ Set a custom font """
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


def lognorm_params(par1, par2):
    """
    Given the mean and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    mean = np.log(par1 ** 2 / np.sqrt(par2 ** 2 + par1 ** 2))  # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2 ** 2 / par1 ** 2 + 1))  # Computes sigma for the underlying normal distribution

    scale = np.exp(mean)
    shape = sigma
    return shape, scale



class dwelltime_by_genotype(hpv.Analyzer):
    '''
    Determine the age at which people with cervical cancer were causally infected and
    time spent between infection and cancer.
    '''

    def __init__(self, start_year=None, **kwargs):
        super().__init__(**kwargs)
        self.start_year = start_year
        self.years = None

    def initialize(self, sim):
        super().initialize(sim)
        self.years = sim.yearvec
        if self.start_year is None:
            self.start_year = sim['start']
        self.age_causal = []
        self.age_cancer = []
        self.dwelltime = dict()
        for state in ['precin', 'cin1', 'cin2', 'cin3', 'total']:
            self.dwelltime[state] = dict()
            for gtype in range(sim['n_genotypes']):
                self.dwelltime[state][gtype] = []

    def apply(self, sim):
        if sim.yearvec[sim.t] >= self.start_year:
            cancer_genotypes, cancer_inds = (sim.people.date_cancerous == sim.t).nonzero()
            if len(cancer_inds):
                current_age = sim.people.age[cancer_inds]
                date_exposed = sim.people.date_exposed[cancer_genotypes, cancer_inds]
                date_cin1 = sim.people.date_cin1[cancer_genotypes, cancer_inds]
                date_cin2 = sim.people.date_cin2[cancer_genotypes, cancer_inds]
                date_cin3 = sim.people.date_cin3[cancer_genotypes, cancer_inds]
                hpv_time = (date_cin1 - date_exposed) * sim['dt']
                cin1_time = (date_cin2 - date_cin1) * sim['dt']
                cin2_time = (date_cin3 - date_cin2) * sim['dt']
                cin3_time = (sim.t - date_cin3) * sim['dt']
                total_time = (sim.t - date_exposed) * sim['dt']
                self.age_causal += (current_age - total_time).tolist()
                self.age_cancer += current_age.tolist()
                for gtype in range(sim['n_genotypes']):
                    gtype_inds = hpv.true(cancer_genotypes == gtype)
                    self.dwelltime['precin'][gtype] += hpv_time[gtype_inds].tolist()
                    self.dwelltime['cin1'][gtype] += cin1_time[gtype_inds].tolist()
                    self.dwelltime['cin2'][gtype] += cin2_time[gtype_inds].tolist()
                    self.dwelltime['cin3'][gtype] += cin3_time[gtype_inds].tolist()
                    self.dwelltime['total'][gtype] += total_time[gtype_inds].tolist()
        return

    def finalize(self, sim=None):
        super().initialize(sim)

