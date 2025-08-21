# BHARP_code

This repository contains code for the simulations in the manuscript:  
*"Identifying Treatment Effect Heterogeneity with Bayesian Hierarchical Adjustable Random Partition (BHARP) in Adaptive Enrichment Trials"*.

Maintained by Xianglin Zhao.

## Structure

- `functions/`  
  Core functions and model implementations:
  - `BHARP.cpp` — rjMCMC implementation of BHARP 
  - `BHARP_model.R` — R wrapper for BHARP  
  - `BHM_model.R`, `IND_model.R`, `BLAST_model.R`, `BART_model.R` — comparator models  
  - `helper_metric.R`, `helper_trial.R` — helper functions

- `scripts/`  
  Scripts for running simulation studies.  
  - `BHARP_sim.R` generates 500 datasets for each scenario and analyzes them with BHARP.  
  - Other method scripts (e.g., BHM, IND, BLAST, BART) do not generate datasets; they analyze the same datasets created by `BHARP_sim.R`.


- `results/`  
  Output directory, stores files written by samplers.  
  Generated datasets are removed but can be reproduced from the provided code.
  Summary results (RMSE, MAE, variance) for each scenario (S1–S9) are provided in `.txt` files (e.g., `S1_BHARP_sim_results.txt`, `S1_BART_sim_results.txt`, etc.).

## Requirements

- **R** (>= 4.0)  
- R packages: `Rcpp`, `RcppArmadillo`, plus others listed in the individual scripts  
- C++14 compiler (for `BHARP.cpp`)

## Usage

1. Run `scripts/BHARP_sim.R` to generate data and run BHARP simulation.  
   This script must be run first, as it creates the simulation datasets.  

2. Run other methods in `scripts/` to process the output in `results/simulation/`.
Each script automatically loads the required functions from `functions/`.
