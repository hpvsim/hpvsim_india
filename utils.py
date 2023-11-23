"""
Utilities
"""

# Imports
import sciris as sc
import numpy as np
import hpvsim as hpv
from scipy.stats import norm, lognorm
import pandas as pd

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
        self.age_causal = dict()
        self.age_cancer = dict()
        self.dwelltime = dict()
        self.median_age_causal = dict()
        for gtype in range(sim['n_genotypes']):
            self.age_causal[gtype] = []
            self.age_cancer[gtype] = []
        for state in ['precin', 'cin', 'total']:
            self.dwelltime[state] = dict()
            for gtype in range(sim['n_genotypes']):
                self.dwelltime[state][gtype] = []

    def apply(self, sim):
        if sim.yearvec[sim.t] >= self.start_year:
            cancer_genotypes, cancer_inds = (sim.people.date_cancerous == sim.t).nonzero()
            if len(cancer_inds):
                current_age = sim.people.age[cancer_inds]
                date_exposed = sim.people.date_exposed[cancer_genotypes, cancer_inds]
                dur_precin = sim.people.dur_precin[cancer_genotypes, cancer_inds]
                dur_cin = sim.people.dur_cin[cancer_genotypes, cancer_inds]
                total_time = (sim.t - date_exposed) * sim['dt']
                for gtype in range(sim['n_genotypes']):
                    gtype_inds = hpv.true(cancer_genotypes == gtype)
                    self.dwelltime['precin'][gtype] += dur_precin[gtype_inds].tolist()
                    self.dwelltime['cin'][gtype] += dur_cin[gtype_inds].tolist()
                    self.dwelltime['total'][gtype] += total_time[gtype_inds].tolist()
                    self.age_causal[gtype] += (current_age[gtype_inds] - total_time[gtype_inds]).tolist()
                    self.age_cancer[gtype] += (current_age[gtype_inds]).tolist()
        return

    def finalize(self, sim=None):
        ''' Convert things to arrays '''
        for gtype in range(sim['n_genotypes']):
            self.median_age_causal[gtype] = np.quantile(self.age_causal[gtype], 0.5)


class age_causal(hpv.Analyzer):
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
        self.age_cin = []


    def apply(self, sim):
        if sim.yearvec[sim.t] >= self.start_year:
            cancer_genotypes, cancer_inds = (sim.people.date_cancerous == sim.t).nonzero()
            if len(cancer_inds):
                current_age = sim.people.age[cancer_inds]
                date_exposed = sim.people.date_exposed[cancer_genotypes, cancer_inds]
                dur_cin = sim.people.dur_cin[cancer_genotypes, cancer_inds]
                total_time = (sim.t - date_exposed) * sim['dt']
                self.age_causal += (current_age - total_time).tolist()
                self.age_cin += (current_age - dur_cin).tolist()
                self.age_cancer += (current_age).tolist()
        return

    def finalize(self, sim=None):
        ''' Convert things to arrays '''
        return


def percentiles_to_pars(x1, p1, x2, p2):
    """ Find the parameters of a normal distribution where:
            P(X < p1) = x1
            P(X < p2) = x2
    """
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)

    location = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = (x2 - x1) / (p2ppf - p1ppf)
    return location, scale


def logn_percentiles_to_pars(x1, p1, x2, p2):
    """ Find the parameters of a lognormal distribution where:
            P(X < p1) = x1
            P(X < p2) = x2
    """
    x1 = np.log(x1)
    x2 = np.log(x2)
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)
    s = (x2 - x1) / (p2ppf - p1ppf)
    mean = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = np.exp(mean)
    return s, scale


def read_debut_data(dist_type='lognormal'):
    """
    Read in dataframes taken from DHS and return them in a plot-friendly format,
    optionally saving the distribution parameters
    """
    df1 = pd.read_csv('data/afs_dist.csv')
    df2 = pd.read_csv('data/afs_median.csv')

    # Deal with median data
    df2['y'] = 50

    # Rearrange data into a plot-friendly format
    dff = {}
    rvs = {'Women': {}, 'Men': {}}

    for sex in ['Women', 'Men']:

        dfw = df1[['Country', f'{sex} 15', f'{sex} 18', f'{sex} 20', f'{sex} 22', f'{sex} 25', f'{sex} never']]
        dfw = dfw.melt(id_vars='Country', value_name='Percentage', var_name='AgeStr')

        # Add values for proportion ever having sex
        countries = dfw.Country.unique()
        n_countries = len(countries)
        vals = []
        for country in countries:
            val = 100-dfw.loc[(dfw['AgeStr'] == f'{sex} never') & (dfw['Country'] == country) , 'Percentage'].iloc[0]
            vals.append(val)

        data_cat = {'Country': countries, 'AgeStr': [f'{sex} 60']*n_countries}
        data_cat["Percentage"] = vals
        df_cat = pd.DataFrame.from_dict(data_cat)
        dfw = pd.concat([dfw,df_cat])

        conditions = [
            (dfw['AgeStr'] == f"{sex} 15"),
            (dfw['AgeStr'] == f"{sex} 18"),
            (dfw['AgeStr'] == f"{sex} 20"),
            (dfw['AgeStr'] == f"{sex} 22"),
            (dfw['AgeStr'] == f"{sex} 25"),
            (dfw['AgeStr'] == f"{sex} 60"),
        ]
        values = [15, 18, 20, 22, 25, 60]
        dfw['Age'] = np.select(conditions, values)

        dff[sex] = dfw

        res = dict()
        res["location"] = []
        res["par1"] = []
        res["par2"] = []
        res["dist"] = []
        for pn,country in enumerate(countries):
            dfplot = dfw.loc[(dfw["Country"] == country) & (dfw["AgeStr"] != f'{sex} never') & (dfw["AgeStr"] != f'{sex} 60')]
            x1 = 15
            p1 = dfplot.loc[dfplot["Age"] == x1, 'Percentage'].iloc[0] / 100
            x2 = df2.loc[df2["Country"]==country,f"{sex} median"].iloc[0]
            p2 = .50
            res["location"].append(country)
            res["dist"].append(dist_type)

            if dist_type=='normal':
                loc, scale = percentiles_to_pars(x1, p1, x2, p2)
                rv = norm(loc=loc, scale=scale)
                res["par1"].append(loc)
                res["par2"].append(scale)
            elif dist_type=='lognormal':
                s, scale = logn_percentiles_to_pars(x1, p1, x2, p2)
                rv = lognorm(s=s, scale=scale)
                res["par1"].append(rv.mean())
                res["par2"].append(rv.std())

            rvs[sex][country] = rv

        pd.DataFrame.from_dict(res).to_csv(f'data/sb_pars_{sex.lower()}_{dist_type}.csv')

    return countries, dff, df2, rvs


class AFS(hpv.Analyzer):
    def __init__(self, bins=None, cohort_starts=None, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins or np.arange(12,31,1)
        self.cohort_starts = cohort_starts
        self.binspan = self.bins[-1]-self.bins[0]

    def initialize(self, sim):
        super().initialize()
        if self.cohort_starts is None:
            first_cohort = sim['start'] + sim['burnin'] - 5
            last_cohort = sim['end']-self.binspan
            self.cohort_starts = sc.inclusiverange(first_cohort, last_cohort)
            self.cohort_ends = self.cohort_starts+self.binspan
            self.n_cohorts = len(self.cohort_starts)
            self.cohort_years = np.array([sc.inclusiverange(i,i+self.binspan) for i in self.cohort_starts])

        self.prop_active_f = np.zeros((self.n_cohorts,self.binspan+1))
        self.prop_active_m = np.zeros((self.n_cohorts,self.binspan+1))

    def apply(self, sim):
        if sim.yearvec[sim.t] in self.cohort_years:
            cohort_inds, bin_inds = sc.findinds(self.cohort_years, sim.yearvec[sim.t])
            for ci,cohort_ind in enumerate(cohort_inds):
                bin_ind = bin_inds[ci]
                bin = self.bins[bin_ind]

                conditions_f = sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin) * sim.people.level0
                denom_inds_f = hpv.true(conditions_f)
                num_conditions_f = conditions_f * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_f = hpv.true(num_conditions_f)
                self.prop_active_f[cohort_ind,bin_ind] = len(num_inds_f)/len(denom_inds_f)

                conditions_m = ~sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin)
                denom_inds_m = hpv.true(conditions_m)
                num_conditions_m = conditions_m * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_m = hpv.true(num_conditions_m)
                self.prop_active_m[ci,bin_ind] = len(num_inds_m)/len(denom_inds_m)
        return


class prop_married(hpv.Analyzer):
    def __init__(self, bins=None, years=None, includelast=True, yearstride=5, binspan=5, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins or np.arange(15, 50, binspan)
        self.years = years
        self.dfs = sc.autolist()
        self.df = None
        self.includelast = includelast
        self.yearstride = yearstride
        self.binspan = binspan

    def initialize(self, sim):
        super().initialize()
        if self.years is None:
            start = sim['start'] + sim['burnin']
            end = sim['end']
            self.years = np.arange(start, end, self.yearstride)
            if self.includelast:
                if end not in self.years:
                    self.years = np.append(self.years, end)

    def apply(self, sim):
        if sim.yearvec[sim.t] in self.years:

            conditions = dict()
            for ab in self.bins:
                conditions[ab] = (sim.people.age >= ab) & (sim.people.age < ab+self.binspan) & sim.people.alive & sim.people.is_female & sim.people.level0

            prop_married = sc.autolist()
            for age_cond in conditions.values():
                num_condition = age_cond & (sim.people.current_partners[0,:]>0)
                prop_married += len(hpv.true(num_condition))/len(hpv.true(age_cond))

            d = dict(age=self.bins, val=prop_married)
            df = pd.DataFrame().from_dict(d)
            df['year'] = sim.yearvec[sim.t]
            self.dfs += df

    def finalize(self, sim):
        self.df = pd.concat(self.dfs)


class cum_dist(hpv.Analyzer):
    """
    Jamie's analyzer for determining distribution of time to clearance, persistence, pre-cancer and cancer
    Debugging differences
        - subsetting men?
    """

    def __init__(self, start_year=None, **kwargs):
        super().__init__(**kwargs)
        self.start_year = start_year

    def initialize(self, sim):
        super().initialize(sim)
        if self.start_year is None:
            self.start_year = sim['start']
        self.dur_to_clearance = []
        self.dur_to_cin = []
        self.dur_to_cancer = []
        self.total_infections = 0

    def apply(self, sim):
        if sim.yearvec[sim.t] >= self.start_year:
            # inf_genotypes, inf_inds = (sim.people.date_exposed == sim.t).nonzero()
            inf_genotypes, inf_inds = ((sim.people.date_exposed == sim.t) & (sim.people.sex==0)).nonzero()
            self.total_infections += len(inf_inds)
            if len(inf_inds):
                infs_that_progress_bools = hpv.utils.defined(sim.people.date_cin[inf_genotypes, inf_inds])
                infs_that_progress_inds = hpv.utils.idefined(sim.people.date_cin[inf_genotypes, inf_inds], inf_inds)
                infs_to_cancer_bools = hpv.utils.defined(sim.people.date_cancerous[inf_genotypes, inf_inds])
                infs_to_cancer_inds = hpv.utils.idefined(sim.people.date_cancerous[inf_genotypes, inf_inds], inf_inds)
                infs_that_clear_bools = hpv.utils.defined(sim.people.date_clearance[inf_genotypes, inf_inds])
                infs_that_clear_inds = hpv.utils.idefined(sim.people.date_clearance[inf_genotypes, inf_inds], inf_inds)

                dur_to_clearance = (sim.people.date_clearance[inf_genotypes[infs_that_clear_bools], infs_that_clear_inds] - sim.t)*sim['dt']
                dur_to_cin = (sim.people.date_cin[inf_genotypes[infs_that_progress_bools], infs_that_progress_inds] - sim.t)*sim['dt']
                dur_to_cancer = (sim.people.date_cancerous[inf_genotypes[infs_to_cancer_bools], infs_to_cancer_inds] - sim.t)*sim['dt']

                self.dur_to_clearance += dur_to_clearance.tolist()
                self.dur_to_cin += dur_to_cin.tolist()
                self.dur_to_cancer += dur_to_cancer.tolist()


class outcomes_by_year(hpv.Analyzer):
    def __init__(self, start_year=1960, **kwargs):
        super().__init__(**kwargs)
        self.start_year = start_year
        self.interval = 1
        self.durations = np.arange(0, 51, self.interval)
        result_keys = ['cleared', 'persisted', 'progressed', 'cancer', 'dead', 'total']
        self.results = {rkey: np.zeros_like(self.durations) for rkey in result_keys}

    def initialize(self, sim):
        super().initialize(sim)
        if self.start_year is None:
            self.start_year = sim['start']

    def apply(self, sim):
        if sim.yearvec[sim.t] == self.start_year:

            conds = (sim.people.date_exposed == sim.t) & (sim.people.sex == 0)  # Get people exposed on this step
            scale = sim.people.scale

            for idd, dd in enumerate(self.durations):

                dur_overall = sim.people.dur_infection + sim.people.dur_cancer

                # cleared = (conds
                #            & ((sim.people.dur_infection <= dd) | (sim.people.date_clearance == sim.t+dd))  # Infections shorted than dd
                #            & (~np.isnan(sim.people.date_clearance))  # Have a clearance date i.e. no cancer
                #            & (np.isnan(sim.people.date_cin) | (sim.people.dur_precin > dd)))  # either no CIN or CIN later
                cleared = (conds
                           & (sim.people.date_clearance <= (sim.t+dd*sim['dt']))  # Infections shorted than dd
                           & (np.isnan(sim.people.date_cin) | (sim.people.dur_precin > dd)))  # either no CIN or CIN later
                persisted = (conds & (sim.people.dur_infection > dd)
                           & (sim.people.dur_precin > dd))
                progressed = (conds & (sim.people.dur_infection > dd)
                              & ~np.isnan(sim.people.date_cin)
                              & (sim.people.dur_precin <= dd)
                              & (np.isnan(sim.people.date_cancerous) | (sim.people.dur_infection > dd)))
                cancer = (conds & ~np.isnan(sim.people.date_cancerous)
                          & (sim.people.dur_infection <= dd)
                          & (dur_overall > dd))
                dead = conds & (dur_overall <= dd)

                # if dd == 30:
                #     import traceback;
                #     traceback.print_exc();
                #     import pdb;
                #     pdb.set_trace()

                cleared_inds = hpv.true(cleared)
                persisted_inds = hpv.true(persisted)
                progressed_inds = hpv.true(progressed)
                cancer_inds = hpv.true(cancer)
                dead_inds = hpv.true(dead)
                derived_total = len(cleared_inds) + len(persisted_inds) + len(progressed_inds) + len(cancer_inds) + len(dead_inds)

                if derived_total != len(hpv.true(conds)):
                    import traceback;
                    traceback.print_exc();
                    import pdb;
                    pdb.set_trace()
                    errormsg = "Something is wrong!"
                    raise ValueError(errormsg)
                scaled_total = scale[hpv.true(conds)].sum()
                self.results['cleared'][idd] += scale[cleared_inds].sum()
                self.results['persisted'][idd] += scale[persisted_inds].sum()
                self.results['progressed'][idd] += scale[progressed_inds].sum()
                self.results['cancer'][idd] += scale[cancer_inds].sum()
                self.results['dead'][idd] += scale[dead_inds].sum()
                self.results['total'][idd] += scaled_total


# %% Run as a script
if __name__ == '__main__':

    countries, dff, df2, rvs = read_debut_data(dist_type='lognormal')

