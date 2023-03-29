# ABM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7271552.svg)](https://doi.org/10.5281/zenodo.7271552)

This readme file contains an overview of the folders with the datasets and programs used to run the models in the article "Economic forecasting with an agent-based model" (Poledna et al. 2023). The programs are written in MATLAB, and data is stored as binary MATLAB files ( `-MAT` files). In the spirit of Dynare, the model is implemented almost as it is described in the manuscript. Thus, the Online Appendix of the article with the full description of the model also documents the implementation. In particular, the function abm.m implements the model with an almost one-to-one correspondence with the formal description and the syntax in the Online Appendix A.

## Basic Simulation

To run a basic simulation with the model in MATLAB:

1. Go to the folder `./model/` in MATLAB.
2. Run the command `run_abm`, this will run out-of-sample forecasts from 2010Q1 to 2019Q4. Or
3. Run the command `[xxx, xxx, xxx, ...] = simulateabm(year, quarter, seed, predictors)`

[nominalgdp,realgdp]=simulateabm(2010,4,1,0) simulates the ABM from 2010:Q4 with seed 1 for the random number generator and returns twelve quarters ahead out-of-sample forecasts for nominal GDP and real GDP and other variables (please refer to simulateabm.m for more information).

## Description of Main Folders

### Calibration

Contains scripts and functions to calibrate model parameters. The script `set_parameters_and_initial_conditions.m` loads the data and calibrates parameters and initial conditions to 40 different reference quarters from the first quarter of 2010 to the third quarter of 2019. These estimates are saved in the `./model/` folder under `parameters` and `initial_conditions` folders, respectively.

### Data

Contains all the required time series data stored as `-MAT` files as well as simulation results from all models used for calibration, estimation, and validation of the models.

### Model

Contains scripts and functions of the ABM (scaled version 1:10000). The function `abm.m` implements the model close to the formal description. The functions `simulateabm.m` and simulate `simulateabmmc.m` run individual and Monte Carlo simulations of the ABM respectively. The script `runabm.m` runs Monte Carlo simulations for the 40 reference quarters from `2010:Q1` to `2019:Q4`. The scripts `import_abm.m` and `import_abmx.m` save ABM forecasts to `./data/abm/` and `./data/abmx/` respectively, with `abmx` conditioning results on some variables being exogenous.

## References

Poledna, S., Miess, M. G., Hommes, C., & Rabitsch, K. (2023). Economic forecasting with an agent-based model. European Economic Review, 151, 104306.
