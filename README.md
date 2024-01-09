# HPVsim model for India

This repository includes the code for creating and calibrating a model of HPV transmission and progression in India. The repository is organized to support reproducing the figures and analyses of the HPVsim methods manuscript.


## Installation

If HPVsim is already installed (`pip install hpvsim`), the only other required dependency is `seaborn`.


## Organization

The repository is organized as follows:

### Inputs
- `data` includes all the input datafiles used for making and validating the model.
- `behavior_inputs.py` contains additional behavioral input parameters.

### Running scripts
- `run_sim.py` contains separate methods for running the sim, extracting sexual behavior, running and plotting a calibration, and running multiple parsets. In order to use the plotting scripts described below, these functions must be run first.

### Plotting scripts

There are separate scripts for plotting each figure in the HPVsim methods manuscript. Specifically:

#### `plot_fig2.py`
 - This script can be used to reproduce Figure 2.

#### `plot_figS1.py` 
- This script can be used to reproduce Figure S1.

#### `plot_figS2.py` 
- This script can be used to reproduce Figure S2.

#### `plot_figS3.py` 
- This script can be used to reproduce Figure S3.

### Additional utility scripts
- `utils.py` contains utilities for numerical calculations, processing data, and creating plots.
- `plot_degree.py` plots the network degree distribution for casual partnerships.
