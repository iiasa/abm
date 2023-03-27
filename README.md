# Readme File

This readme file contains an overview of the folders with the datasets and program used to run the models. The programs are written in MATLAB, and data is stored as binary MATLAB files ( `-MAT` files). In the spirit of Dynare, the model is implemented almost as it is described in the manuscript. Thus, the Online Appendix with the full description of the model also documents the implementation. In particular, the function abm.m implements the model with an almost one-to-one correspondence with the formal description and the syntax in the Online Appendix A.

## Basic Simulation

To run a basic simulation with the model in MATLAB:

1. Go to the folder `./model/` in MATLAB.
2. Run the command `run_abm`, this will run out-of-sample forecasts from 2010Q1 to 2019Q4. Or
3. Run the commnad `[xxx, xxx, xxx] = run(year, quarter, seed ,predictors)`

[nominal\_gdp,real\_gdp,xxx variables]=simulate\_abm(2010,4,1,0)} simulates the ABM from 2010:Q4 with seed 1 for the random number generator and returns twelve quarters ahead out-of-sample forecasts for nominal GDP and real GDP and other variables (please refer to simulate\_abm.m for more information).

## Description of Main Folders

### Calibration

Contains scripts and functions to calibrate model parameters. The script `set_parameters_and_initial_conditions.m` loads the data and calibrates parameters and initial conditions to 40 different reference quarters from the first quarter of 2010 to the third quarter of 2019. These estimates are saved in the `./model/` folder under `parameters` and `initial_conditions` folders, respectively.

### Data

Contains all the required time series data stored as `-MAT` files as well as simulation results from all models used for calibration, estimation, and validation of the models.

### Model

Contains scripts and functions of the ABM. The functions `abm.m` defines an individual simulation of the ABM whereas the script `run.m` runs one iteration of the model for the 40 reference quarters from `2010:Q1` to `2019:Q4`. The scripts `import_abm.m` and `import_abmx.m` saves ABM forecasts to `./data/abm/` and `./data/abmx/` respectively, with `abmx` conditioning results on some variables being exogenous.
