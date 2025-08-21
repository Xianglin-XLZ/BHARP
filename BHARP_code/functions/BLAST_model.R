# File: BLAST_model.R
# Purpose:
#   This script defines the BLAST model used as a baseline comparator in the simulation
#   of treatment effect heterogeneity. BLAST is a fixed finite mixture model (FMM) model.
#   The number of components (q) is selected via Deviance Information Criterion (DIC).
#
# Contents:
#   - BLAST_code: Stan model specification for the FMM with q components.
#   - BLAST_stan_model: Compiled Stan model using `rstan::stan_model()`.
#   - BLAST_analysis(): Fit the BLAST model with q = 1, 2, 3 and select best q via DIC.
#   - dic(): Compute the DIC for model comparison based on posterior log-likelihood.
#
# Dependencies:
#   - Requires the `rstan` and `posterior` packages.
#   - Assumes data input is in the form of a list with outcome vector Y and subgroup index k.


# ----BLAST----

BLAST_code<-"

data {
  int<lower=0> N;                 // total sample size 
  int<lower=1> Nsubgrp;           // number of subgroups 
  vector[N] Y;                    // outcomes
  int<lower=1,upper=Nsubgrp> kIndex[N]; // subgroup index for each observation
  
  int<lower=1,upper=Nsubgrp> q; //number of mixing component
  
  //hyperparam  for FMM variances
  real a_between;
  real b_between; 
  real a_within;
  real b_within;
  
  // hyperparam for cell precision
  real a_cell;    // shape
  real b_cell;    // rate
  
  //hyperparam for beta
  real c_beta;
  real p_beta; 

}

parameters {
  real beta;
  real<lower=0> cpre;
  real<lower=0> tau;
  simplex[q] w;
  vector<lower=0>[q] sigma;
  vector[q] mu;
  
  vector[Nsubgrp] theta;
  
  
}

model{
  //prior 
   beta ~ normal(c_beta, inv_sqrt(p_beta));
   cpre ~ gamma(a_cell, b_cell);
   tau ~ gamma(a_between,b_between);
   mu ~ normal(beta, inv_sqrt(tau));
   sigma ~inv_gamma(a_within,b_within);
   w ~ dirichlet(rep_vector(1.0, q));
   
  //mixture
  for (k in 1:Nsubgrp) {
    vector[q] lp_mix;
    for (t in 1:q)
      lp_mix[t] = log(w[t]) + normal_lpdf(theta[k] | mu[t], sqrt(sigma[t]));
    target += log_sum_exp(lp_mix);  
  }
  
  Y ~ normal(theta[kIndex], inv_sqrt(cpre));
}

"

BLAST_stan_model <- rstan::stan_model(model_code = BLAST_code)
BLAST_analysis<-function(Data,n_subgrp, hyperparam){
  
  standata <- list(
    N      = length(Data$Y),
    Nsubgrp= n_subgrp,
    Y      = Data$Y,
    kIndex = Data$k,
    
    a_cell     = hyperparam$a_cell,
    b_cell     = hyperparam$b_cell,
    a_within =hyperparam$a_within,
    b_within =hyperparam$b_within,
    a_between =hyperparam$a_between,
    b_between =hyperparam$b_between,
    c_beta=hyperparam$c_beta,
    p_beta=hyperparam$p_beta
  )
  
  fit_list <- lapply(1:3, function(q_now) {
    rstan::sampling(BLAST_stan_model,             
                    data = c(standata, list(q = q_now)),
                    chains = 4, iter = 2000, seed = 2025)
  })
  
  
  dic_tbl <- do.call(rbind,
                     Map(function(f, q_now) cbind(q = q_now, dic(f)),
                         fit_list, 1:3))
  best_q <- dic_tbl$q[ which.min(dic_tbl$DIC) ]
  
  samples <- rstan::extract(fit_list[[best_q]])
  return(
    list(q=best_q,theta=samples$theta)
    
  )
  
}




dic <- function(fit) {
  log_lik <- as_draws_matrix(fit, variable = "log_lik")
  D_bar <- -2 * mean(rowSums(log_lik))
  D_hat <- -2 * sum(colMeans(log_lik))
  DIC <- 2 * D_bar - D_hat
  data.frame(D_bar, D_hat, DIC)
}



