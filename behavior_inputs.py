"""
Compilation of sexual behavior data and assumptions
"""


#%% Initialization

import numpy as np

#%% LAYER PROBS
layer_probs = dict(
    m=np.array([
        [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75],
        [0, 0, 0.05, 0.25, 0.60, 0.80, 0.95, 0.80, 0.80, 0.65, 0.55, 0.40, 0.40, 0.40, 0.40, 0.40],  # Share f
        [0, 0, 0.01, 0.05, 0.10, 0.70, 0.90, 0.90, 0.90, 0.90, 0.80, 0.60, 0.60, 0.60, 0.60, 0.60]]  # Share m
    ),
    c=np.array([
        [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75],
        [0, 0, 0.10, 0.50, 0.60, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.01, 0.01],  # Share f
        [0, 0, 0.10, 0.20, 0.25, 0.35, 0.40, 0.70, 0.90, 0.90, 0.95, 0.95, 0.70, 0.30, 0.10, 0.10]],  # Share m
    ),
)

#%% PARTNERS
m_partners = dict(
    m=dict(dist='poisson1', par1=0.01),
    c=dict(dist='poisson1', par1=0.1),
)
f_partners = dict(
    m=dict(dist='poisson1', par1=0.01),
    c=dict(dist='neg_binomial', par1=2, par2=0.025),
)

