# BHARP

Code for **Bayesian Hierarchical Adjustable Random Partition (BHARP)**, a Bayesian method for borrowing information across subgroups while allowing the subgroup partition to adapt to the observed data.

This repository contains the BHARP sampler, simulation studies, and adaptive trial design evaluations. It also includes implementations of four comparison approaches:

- BLAST
- Bayesian hierarchical model (BHM)
- Independent subgroup analysis (IND)
- Bayesian additive regression trees (BART)

## Repository structure

```text
BHARP_release/
├── functions/   # Model implementations and helper functions
├── scripts/     # Simulation and adaptive trial design scripts
└── results/     # Output directory created/used by the scripts
```

The BHARP MCMC and reversible-jump MCMC sampler is implemented in `functions/BHARP_theta.cpp` and compiled from R with `Rcpp::sourceCpp()`.

## Requirements

The code requires R and a C++14-compatible compiler. Install the R packages used by the main and comparison methods:

```r
install.packages(c(
  "Rcpp",
  "RcppArmadillo",
  "dplyr",
  "purrr",
  "tibble",
  "rstan",
  "BART"
))
```

The `parallel` package is included with R.

> **Platform note:** The simulation scripts use FORK-based parallel processing and are intended for macOS and Linux. They will require modification to run on Windows.

## Usage

Download this repository and update the following line near the top of each script you intend to run:

```r
project_root <- "YOURPATH/BHARP_release"
```

Replace the placeholder with the absolute path to the downloaded repository.

### Estimation simulation study

Simulation settings, including the number of replications, scenarios, sample sizes, and hyperparameters, are defined in:

```text
scripts/simulation_config.R
```

Run the BHARP simulation first:

```r
source("scripts/BHARP_sim.R")
```

This generates the simulated datasets and fits BHARP. Results are written under:

```text
results/simulation/<scenario_name>/
```

The comparison methods use the same generated datasets and can then be run separately:

```r
source("scripts/BLAST_sim.R")
source("scripts/BHM_sim.R")
source("scripts/IND_sim.R")
source("scripts/BART_sim.R")
```

### Adaptive trial design evaluation

Trial settings, including the true cell means, interim sample sizes, decision boundaries, probability thresholds, and hyperparameters, are defined in:

```text
scripts/trial_design_config.R
```

Run the trial simulation for the desired method:

```r
source("scripts/BHARP_trialsim.R")
source("scripts/BLAST_trialsim.R")
source("scripts/BHM_trialsim.R")
source("scripts/IND_trialsim.R")
source("scripts/BART_trialsim.R")
```

Outputs are written under:

```text
results/trial_design/
```

## Reproducibility notes

- The supplied scripts use fixed random seeds.
- The default configurations run 500 simulation replications.
- Model output and generated datasets can require considerable disk space.
- Review the configuration files before running the full simulations.

## Citation

If you use this code, please cite the associated BHARP manuscript. 

```bibtex
@article{zhao2025identifying,
  title={Identifying treatment effect heterogeneity with Bayesian hierarchical adjustable random partition in adaptive enrichment trials},
  author={Zhao, Xianglin and Golchi, Shirin and Gouin, Jean-Philippe and Dasgupta, Kaberi},
  journal={arXiv preprint arXiv:2508.16523},
  year={2025}
}
```
## License


This project is licensed under the MIT License. See `LICENSE` for details.
