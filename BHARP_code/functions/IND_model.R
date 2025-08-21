# File: ind_model.R
# Purpose:
#   This script defines the Bayesian Independent (IND) model used for comparison
#   in simulation studies. The IND model assumes independent treatment effects
#   across all arm-subgroup combinations without information sharing.
#
# Contents:
#   - IND_code:    Stan model specification for the independent model.
#   - IND_stan_model: Compiled Stan model using `rstan::stan_model()`.
#   - IND_analysis():  Fit the IND model to a given dataset and return posterior samples.
#   - IND_simulate_one_trial(): Simulate an adaptive multi-stage clinical trial
#                               using the IND model and predefined decision rules.
#
# Dependencies:
#   - Requires the `rstan` package to compile and run the Stan model.
#   - Depends on external utility functions such as:
#       - recruit_patients()         (for stage-wise patient generation)
#       - deactivate_decision(), final_decision() (for adaptive decision logic)



# ----Bayesian independent----
IND_code <- "
  data {
    int<lower=0> N;            // total sample size
    int<lower=1> Narm;         // arm number
    int<lower=1> Nsubgrp;         // subgroup number
    vector[N] Y;               // outcome 
    int<lower=1,upper=Narm> iIndex[N];             //  arm index i
    int<lower=1,upper=Nsubgrp> kIndex[N];          // subgroup index k

    real c_thetaind;                  // theta prior mean
    real p_thetaind;                   // theta prior precision
    real a_cell;               // tau prior shape
    real b_cell;               // tau  prior rate
  }
  parameters {
    matrix[Narm, Nsubgrp] theta;     // cell center
    real<lower=0> cpre;           //  shared cell precision
  }
 
  model {
    // prior
    for (i in 1:Narm) {
      for (k in 1:Nsubgrp) {
        theta[i, k] ~ normal(c_thetaind, inv_sqrt(p_thetaind));
      }
    }
    cpre ~ gamma(a_cell, b_cell);

    for (n in 1:N) {
      Y[n] ~ normal(theta[iIndex[n], kIndex[n]], inv_sqrt(cpre));
    }
  }
  
generated quantities {
  vector[Narm * Nsubgrp] theta_vec;

  for (i in 1:Narm) {
    for (k in 1:Nsubgrp) {
      theta_vec[(i - 1) * Nsubgrp + k] = theta[i, k];  // row-wise flattening
    }
  }
}
 
  "
IND_stan_model <- rstan::stan_model(model_code = IND_code)
IND_analysis<-function(Data,n_subgrp,n_arm, hyperparam){

  standata <- list(
    N      = length(Data$Y),
    Narm   = n_arm,
    Nsubgrp= n_subgrp,
    Y      = Data$Y,
    iIndex = Data$i,
    kIndex = Data$k,
    
    c_thetaind = hyperparam$c_thetaind,
    p_thetaind = hyperparam$p_thetaind,
    a_cell     = hyperparam$a_cell,
    b_cell     = hyperparam$b_cell
  )

  fit <- rstan::sampling(
    IND_stan_model,
    data   = standata,
    iter   = 2500, chains = 4, warmup=500
  )
  samples <- rstan::extract(fit)
  return(
    samples$theta_vec
  )
  
}



IND_simulate_one_trial <- function(
    trial_index, N_subgrp, N_arm,
    L,                   # number of analyses
    TotalSampleSize,     # A vector containing total sample size at each analysis
    Theta,               # Real Theta
    hyperparam_list,  
    bound_fut, bound_eff,          # futility / efficacy boundaries
    P_fut, P_eff,          # futility / efficacy thresholds
    base_dir = getwd()     # base directory to save all simulation results
){
  
  # create folder 
  sim_dir <- file.path(base_dir, paste0("Sim", trial_index))
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  }
  results_dir <- file.path(sim_dir, "Results")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  }
  data_dir <- file.path(sim_dir, "Data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # initialize containers
  Status_Log <- matrix(NA, nrow=L, ncol=N_subgrp*N_arm)
  cell_actflag <- matrix(TRUE,  nrow=N_arm, ncol=N_subgrp)
  cell_final_efficacy <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  cell_final_futility <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  Dat <- data.frame(i=integer(0), k=integer(0), Y=numeric(0))
  
  # the vectors to return
  thetamedian <- rep(NA, N_subgrp*N_arm)
  thetaq750   <- rep(NA, N_subgrp*N_arm)
  thetaq250   <- rep(NA, N_subgrp*N_arm)
  final_efficacy_vec <- rep(NA, N_subgrp*N_arm)
  final_futility_vec <- rep(NA, N_subgrp*N_arm)
  
  # analyze by stage
  for (l in seq_len(L)) {
    # recruit and save data
    Dat <- recruit_patients(
      dataset = Dat,
      cell_active = cell_actflag,
      targetN = TotalSampleSize[l],
      n_subgrp = N_subgrp,
      n_arm = N_arm,
      Theta = Theta,
      method = "balanced"
    )
    write.csv(Dat, file.path(data_dir, paste0("Data_stage", l, ".csv")))
    
    
    theta_l <- IND_analysis(Data=Dat,
                            n_subgrp=N_subgrp,
                            n_arm=N_arm,
                            hyperparam=hyperparam_list)
    
    
    # interim  vs final
    if (l < L) {
      # interim
      decision_l <- deactivate_decision(
        theta_posterior = theta_l,
        cell_active = cell_actflag,
        cell_final_efficacy = cell_final_efficacy,
        cell_final_futility = cell_final_futility,
        futility_bound = bound_fut,
        efficacy_bound = bound_eff,
        futility_threshold = P_fut[l],
        efficacy_threshold = P_eff[l]
      )
      cell_actflag        <- decision_l$Cell_Active_Updated
      cell_final_efficacy <- decision_l$Cell_Final_Efficacy_Updated
      cell_final_futility <- decision_l$Cell_Final_Futility_Updated
      Status_Log[l, ]     <- as.vector(t(decision_l$Status))
      
      # terminate trial if no active cell
      if (all(!cell_actflag)) {
        # save interim results as final
        write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
        
        
        thetamedian <- apply(theta_l, 2, median)
        thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
        thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
        
        
        final_efficacy_vec <- as.vector(t(cell_final_efficacy))
        write.csv(cell_final_efficacy, file.path(sim_dir, "final_efficacy.csv"),row.names = FALSE)
        
        final_futility_vec <- as.vector(t(cell_final_futility))
        write.csv(cell_final_futility, file.path(sim_dir, "final_futility.csv"),row.names = FALSE)
        write.csv(Status_Log, file.path(sim_dir, "status_log.csv"),row.names = FALSE)
        
        
        
        thetamedian_mat<-matrix(thetamedian,nrow = N_arm,ncol = N_subgrp,byrow=TRUE)
        
        cell_sample_size=as.matrix(table(Dat$i, Dat$k))
        #calculate nonfutile cells' weighted average estimation
        avg_effect_nonfutile <- apply(matrix(seq_len(N_arm)), 1, function(i) {
          nf_idx <- which(!cell_final_futility[i, ])
          if (length(nf_idx)==0) {
            NA
          } else {
            w <- cell_sample_size[i, nf_idx]
            sum(thetamedian_mat[i, nf_idx] * w) / sum(w)
          }
        })
        
        
        #calculate effective cells' average estimation
        avg_effect_effective <- apply(matrix(seq_len(N_arm)), 1, function(i) {
          eff_idx <- which(cell_final_efficacy[i, ])
          if (length(eff_idx)==0) {
            NA
          } else {
            w <- cell_sample_size[i, eff_idx]
            sum(thetamedian_mat[i, eff_idx] * w) / sum(w)
          }
        })
        
        
        
        break  # terminate for l in 1:L 
      }
      
    } 
    else {
      # l == L (final)
      # save theta
      write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
      thetamedian <- apply(theta_l, 2, median)
      thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
      thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
      
      # draw final conclusion for active cells
      decision_l <- final_decision(
        theta_posterior=theta_l,
        cell_active=cell_actflag,
        cell_final_efficacy=cell_final_efficacy,
        cell_final_futility=cell_final_futility,
        futility_bound=bound_fut,
        efficacy_bound=bound_eff,
        futility_threshold=P_fut[l],
        efficacy_threshold=P_eff[l]
      )
      
      cell_final_efficacy <- decision_l$Cell_Final_Efficacy_Updated
      cell_final_futility <- decision_l$Cell_Final_Futility_Updated
      Status_Log[l, ]     <- as.vector(t(decision_l$Status))
      
      final_efficacy_vec <- as.vector(t(cell_final_efficacy))
      write.csv(cell_final_efficacy, file.path(sim_dir, "final_efficacy.csv"),row.names = FALSE)
      
      final_futility_vec <- as.vector(t(cell_final_futility))
      write.csv(cell_final_futility, file.path(sim_dir, "final_futility.csv"),row.names = FALSE)
      write.csv(Status_Log, file.path(sim_dir, "status_log.csv"),row.names = FALSE)
      
      
      thetamedian_mat<-matrix(thetamedian,nrow = N_arm,ncol = N_subgrp,byrow=TRUE)
      cell_sample_size=as.matrix(table(Dat$i, Dat$k))
      #calculate nonfutile cells' weighted average estimation
      avg_effect_nonfutile <- apply(matrix(seq_len(N_arm)), 1, function(i) {
        nf_idx <- which(!cell_final_futility[i, ])
        if (length(nf_idx)==0) {
          NA
        } else {
          w <- cell_sample_size[i, nf_idx]
          sum(thetamedian_mat[i, nf_idx] * w) / sum(w)
        }
      })
      
      
      #calculate effective cells' average estimation
      avg_effect_effective <- apply(matrix(seq_len(N_arm)), 1, function(i) {
        eff_idx <- which(cell_final_efficacy[i, ])
        if (length(eff_idx)==0) {
          NA
        } else {
          w <- cell_sample_size[i, eff_idx]
          sum(thetamedian_mat[i, eff_idx] * w) / sum(w)
        }
      })
      
      
      
    } # end if (l < L) else
    
    
    
  } # end for l in seq_len(L)
  
  # return information for one simulation
  return(list(
    sample_size=TotalSampleSize[l],
    cell_sample_size=as.vector(t(cell_sample_size)),
    theta=theta_l,
    thetamedian=thetamedian,
    thetaq750=thetaq750,
    thetaq250=thetaq250,
    final_efficacy=final_efficacy_vec,
    final_futility=final_futility_vec,
    avg_effect_nonfutile=avg_effect_nonfutile,
    avg_effect_effective=avg_effect_effective
  ))
}

