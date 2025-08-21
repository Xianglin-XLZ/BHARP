# BHARP_code

This repository contains code for the simulations in the manuscript:  
*"Identifying Treatment Effect Heterogeneity with Bayesian Hierarchical Adjustable Random Partition (BHARP) in Adaptive Enrichment Trials"*.

Maintained by Xianglin Zhao

## Structure

- `functions/`  
  Core functions and model implementations:
  - `BHARP.cpp` — rjMCMC implementation of BHARP (split–merge sampler)  
  - `BHARP_model.R` — R wrapper for BHARP  
  - `BHM_model.R`, `IND_model.R`, `BLAST_model.R`, `BART_model.R` — comparator models  
  - Helper functions: `helper_metric.R`, `helper_trial.R`

- `scripts/`  
  Scripts for running simulation studies.

- `results/`  
  Output directory for simulation results.  
  Subfolder `simulation/` stores CSV files written by samplers under each scenario.  
  (The actual datasets are removed but can be reproduced from the provided code.)

## Requirements

- **R** (>= 4.0)  
- R packages: `Rcpp`, `RcppArmadillo`, plus others listed in the individual scripts  
- C++14 compiler (for `BHARP.cpp`)

## Usage

1. Compile the C++ functions via Rcpp:
   ```r
   Rcpp::sourceCpp("functions/BHARP.cpp")
   ```
   
2. Run simulation scripts in scripts/.

3. Results will be written to results/simulation/ as CSV files.
