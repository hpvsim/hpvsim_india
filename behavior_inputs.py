"""
Compilation of sexual behavior data and assumptions
"""


#%% Initialization

import numpy as np

#%% LAYER PROBS
layer_probs = dict(
    m=np.array([
        # Share of females (row 1) and males (row 2) of each age who are married
        [0, 5,  10,    15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],  # Age bracket
        [0, 0, 0.05, 0.25, 0.60, 0.90, 0.95, 0.80, 0.80, 0.65, 0.55, 0.40, 0.40, 0.40, 0.40, 0.40],  # Share f
        [0, 0, 0.01, 0.05, 0.10, 0.70, 0.90, 0.90, 0.90, 0.90, 0.80, 0.60, 0.60, 0.60, 0.60, 0.60]]  # Share m
    ),
    c=np.array([
        # Share of females (row 1) and males (row 2) of each age having casual relationships
        [0, 5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],  # Age bracket
        [0, 0, 0.50, 0.60, 0.60, 0.60, 0.80, 0.80, 0.80, 0.80, 0.80, 0.50, 0.40, 0.10, 0.01, 0.01],  # Share f
        [0, 0, 0.10, 0.40, 0.50, 0.60, 0.80, 0.90, 0.90, 0.80, 0.80, 0.70, 0.50, 0.30, 0.10, 0.10]],  # Share m
    ),
)


#%% PARTNERS
m_partners = dict(
    m=dict(dist='poisson1', par1=0.001),
    c=dict(dist='poisson1', par1=10),
)
f_partners = dict(
    m=dict(dist='poisson1', par1=0.001),
    c=dict(dist='neg_binomial', par1=1, par2=0.2),
)

#%% MIXING
mixing_all = np.array([
    #       0,  5, 10, 15, 20,  25,  30,  35,  40,  45,  50,  55,  60,  65,  70,  75
    [0,     0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [5,     0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [10,    0,  0,  1,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [15,    0,  0,  1,  1,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [20,    0,  0, .5,  1,  1, .01,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [25,    0,  0,  0, .5,  1,   1, .01,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [30,    0,  0,  0,  0, .5,   1,   1, .01,   0,   0,   0,   0,   0,   0,   0,   0],
    [35,    0,  0,  0,  0, .1,  .5,   1,   1, .01,   0,   0,   0,   0,   0,   0,   0],
    [40,    0,  0,  0,  0,  0,  .1,  .5,   1,   1, .01,   0,   0,   0,   0,   0,   0],
    [45,    0,  0,  0,  0,  0,   0,  .1,  .5,   1,   1, .01,   0,   0,   0,   0,   0],
    [50,    0,  0,  0,  0,  0,   0,   0,  .1,  .5,   1,   1,  .01,  0,   0,   0,   0],
    [55,    0,  0,  0,  0,  0,   0,   0,   0,  .1,  .5,   1,   1, .01,   0,   0,   0],
    [60,    0,  0,  0,  0,  0,   0,   0,   0,   0,  .1,  .5,   1,   1, .01,   0,   0],
    [65,    0,  0,  0,  0,  0,   0,   0,   0,   0,   0,  .1,  .5,   1,   1, .01,   0],
    [70,    0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,  .1,  .5,   1,   1, .01],
    [75,    0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,  .1,  .5,   1,   1],
])

mixing = dict()
for k in ['m', 'c']: mixing[k] = mixing_all

