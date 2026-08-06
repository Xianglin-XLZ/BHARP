# -------------------------------------------------------------------------
# IND Model and Trial Simulation Functions
#
# This file implements the independent Bayesian model used as a comparator in
# the simulation study. Treatment effects are estimated independently across
# treatment arms and biomarker subgroups, without borrowing information across
# cells.
#
# Main functions:
#   - IND_analysis(): fits the IND model to a single dataset.
#   - IND_simulate_one_trial(): simulates and analyzes one adaptive trial.
#
# Dependencies:
#   - rstan
#   - recruit_patients(), deactivate_decision(), and final_decision()
#     from helper_trial.R
# -------------------------------------------------------------------------


IND_code <- "
  data {
    int<lower=0> N;            
    int<lower=1> Narm;         
    int<lower=1> Nsubgrp;         
    vector[N] Y;               
    int<lower=1,upper=Narm> iIndex[N];             
    int<lower=1,upper=Nsubgrp> kIndex[N];          

    real c_thetaind;                  
    real p_thetaind;                   
    real a_cell;               
    real b_cell;               
  }
  parameters {
    matrix[Narm, Nsubgrp] theta;     
    real<lower=0> cpre;           
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

IND_analysis <- function( Data, n_subgrp, n_arm, hyperparam, seed=NULL){

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
    iter   = 3000,
    chains = 4, 
    warmup=1000, 
    seed=seed
  )
  samples <- rstan::extract(fit)
  
  theta_post <- samples$theta_vec
  theta_post <- matrix(theta_post, ncol = n_arm * n_subgrp, byrow = FALSE)
  
  if (ncol(theta_post) != n_arm * n_subgrp) {
    stop("theta_vec wrong ncol")
  }
  
  return(theta_post)
  
}



IND_simulate_one_trial <- function(
    trial_index, 
    N_subgrp, 
    N_arm,
    L,                   
    TotalSampleSize,     
    Theta,               
    hyperparam_list,  
    bound_fut, bound_eff,          # futility / efficacy boundaries
    P_fut, P_eff,          # futility / efficacy thresholds
    base_dir,    
    trial_seed = NULL
){

  sim_dir <- file.path(base_dir, paste0("Sim", trial_index))
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  }
  data_dir <- file.path(sim_dir, "Data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  }
  

  Status_Log <- matrix(NA, nrow=L, ncol=N_subgrp*N_arm)
  cell_actflag <- matrix(TRUE,  nrow=N_arm, ncol=N_subgrp)
  cell_final_efficacy <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  cell_final_futility <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  Dat <- data.frame(i=integer(0), k=integer(0), Y=numeric(0))
  

  thetamedian <- rep(NA, N_subgrp*N_arm)
  thetaq750   <- rep(NA, N_subgrp*N_arm)
  thetaq250   <- rep(NA, N_subgrp*N_arm)
  final_efficacy_vec <- rep(NA, N_subgrp*N_arm)
  final_futility_vec <- rep(NA, N_subgrp*N_arm)
  

  for (l in seq_len(L)) {

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

    
    theta_l <- IND_analysis(
      Data=Dat,
      n_subgrp=N_subgrp,
      n_arm=N_arm,
      hyperparam=hyperparam_list,
      seed=trial_seed+l
    )
    
    theta_l <- as.matrix(theta_l)
    if (is.null(dim(theta_l))) {
      stop("theta_l has no dim at trial ", trial_index, ", stage ", l)
    }

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
      
      if (all(!cell_actflag)) {

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
        
        break  
      }
      
    } 
    else {
      # l == L (final)

      write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
      thetamedian <- apply(theta_l, 2, median)
      thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
      thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
      
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
      
      
      
    } 
    
    
    
  } 
  
  # return information for one simulation
  return(list(
    sample_size=TotalSampleSize[l],
    cell_sample_size=as.vector(t(cell_sample_size)),
    theta=theta_l,
    thetamedian=thetamedian,
    thetaq750=thetaq750,
    thetaq250=thetaq250,
    final_efficacy=final_efficacy_vec,
    final_futility=final_futility_vec
  ))
}

