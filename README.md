# Temporary Public Version of *C. difficile* Agent-Based Model

Agent-based model for _Clostridium difficile_ transmission in hospital.


The epidemiological model is implemented in C++ and packaged in a R library named `abmcdiff`.
 

## Directory structure

* `src`: C++ code of the epidemiological model.
* `Rlibrary`: scripts to build the R library `abmcdiff`
* `simul`: R scripts that run simulations using the compiled library `abmcdiff`. 
* `doc`: model documentation

## Running simulations

### R library 

First, the R library `abmcdiff` that wraps the C++ model must be built. Execute `build_library` in the `Rlibrary` folder. 

### Model setup

The architecture of the rooms for each health care setting simulated must be specified in a CSV file. See `room_hosp.csv` for an example. 

Then, the parameters defining the population modelled in the health care setting must be defined in another CSV file. These includes the number of HCW staff, the stay duration distribution parameters, the admission rate, contact rates between individuals, contamination decay rates, etc. See `prm_hosp.csv` for an example.

If there are more than one health care setting modelleed (e.g., a hospital and a long term care facility), the movement between the facilities can be defined in the form of a matrix of transfer rates (with 0 on the diagonal). See `prm_mvt_hcs.csv` for an example. 

Finally, the simulation parameters -- like number of Monte Carlo iterations, the time horizon and time step -- must be defined in a separate CSV file. See `prm_sim.csv` for an example. 

### Running a simulation

The R library `abmcdiff` must be loaded. 

Once all CSV files have been saved to define the model parameters, they can be imported in a R script and loaded as R objects (typically as lists). 

One of the main R/C++ wrapping function, `abmcdiff_one_simul()`, takes those parameter objects as inputs and run the disease transmission simulations in the health care settings defined. See `run.R` for an example. 


## Simulation analysis

The analysis is performed separately, loading the simulation outputs saved in a RData file (e.g., `simul.RData`). 


## Goal of the study

We want to investigate the cost-effectiveness of a potential vaccine against C. difficile  currently in stage 2. Its main feature is to reduce the risk of _symptomatic_ infection (does not provide full immunity), and provides this partial immunity for a limited period of time (2 years?).

There are two vaccination strategies envisaged:
* Vaccinate some patients before admission. Only applies to non-emergency admissions. 
* Vaccinate patients who had a CDI during their stay. The simulation must keep track of them if/when they are re-admitted.

For the cost effectiveness analysis, the costs must be recorded at the individual level, with a date attached (for discounting). 
The costs associated with CDI are:
* isolation
* treatment (drugs)
* potentially longer hospital stay









