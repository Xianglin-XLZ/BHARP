# -------------------------------------------------------------------------
# BHARP Model and Trial Simulation Functions
#
# This file defines functions for fitting the BHARP model and simulating
# one adaptive trial.
#
# Functions:
#   ReadCombineMCMC():
#     Reads and combines MCMC samples stored in multiple CSV files.
#
#   Bharp_analysis():
#     Runs the BHARP MCMC algorithm and combines posterior samples across
#     chains.

#   simulate_one_trial():
#     Simulates recruitment, interim analyses, early stopping decisions,
#     and the final analysis for one adaptive trial.
#
# Dependencies:
#   Bharp() is defined in functions/BHARP_theta.cpp.
#
#   The following functions are defined in functions/helper_cocluster.R:
#     ArmListCoCluster()
#
#   The following functions are defined in functions/helper_trial.R:
#     recruit_patients(), deactivate_decision(), final_decision()
# -------------------------------------------------------------------------

# to read and combine MCMC samples 
ReadCombineMCMC<- function(prefix,header=FALSE) {
  dir_path <- dirname(prefix)  
  file_prefix <- basename(prefix)  
  if (dir_path == ".") {
    dir_path <- getwd()
  }
  file_pattern <- paste0("^", file_prefix, ".*\\.csv$")
  csv_files <- list.files(path = dir_path, pattern = file_pattern, full.names = TRUE)
  
  if (length(csv_files) > 0L) {
    combined_data <- do.call(rbind, lapply(csv_files, function(file) {
      read.csv(file, header = header)
    }))
    return(combined_data)
  } else {
    warning(paste("No files found in", dir_path, "matching pattern:", file_pattern))
    return(NULL)
  }
}

Bharp_analysis<-function(
    l,
    dataset,
    n_subgrp,
    n_arm,
    hyperparam_list,
    results_dir
){

  
Bharp(
  Data_i = dataset$i, 
  Data_k = dataset$k, 
  Data_Y = dataset$Y, 
  ngrp = n_subgrp, 
  narm = n_arm,
  a_cell = hyperparam_list$a_cell,       
  b_cell = hyperparam_list$b_cell,
  a_between = hyperparam_list$a_between, 
  b_between = hyperparam_list$b_between,
  a_within = hyperparam_list$a_within,   
  b_within = hyperparam_list$b_within,
  results_dir = results_dir,
  nsamp=2000, 
  thn=3, 
  burn=1000,
  dataid=l,
  nchain=4
)
  

  thetaPost<-ReadCombineMCMC( paste0(results_dir,"/theta_",l,"-") )
  izPost<-as.matrix(ReadCombineMCMC(paste0(results_dir,"/iz_",l,"-")))
  iqPost<-as.matrix(ReadCombineMCMC(paste0(results_dir,"/iq_",l,"-")))
  
  return(
    list(
      theta_posterior = thetaPost,
      iz_posterior = izPost,
      iq_posterior = iqPost
    )
  )
  
}  




# ---- Simulate one adaptive trial using BHARP ----

simulate_one_trial <- function(
    trial_index, 
    N_subgrp,
    N_arm,
    L,                   
    TotalSampleSize,     
    Theta,               
    hyperparam_list,  
    bound_fut, bound_eff,          # futility / efficacy boundaries
    P_fut, P_eff,          # futility / efficacy thresholds
    base_dir     
){
  

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
    # Fit the BHARP model to the accumulated data.
    analysis_result_l <- Bharp_analysis(
      l=l,
      dataset=Dat,
      n_subgrp=N_subgrp,
      n_arm=N_arm,
      hyperparam_list =hyperparam_list, 
      results_dir = results_dir  
    )
    
    theta_l <- analysis_result_l$theta_posterior
    iz_l    <- analysis_result_l$iz_posterior
    
    # ---- Interim analysis ----
    if (l < L) {

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
      
      # Terminate the trial when no cells remain active.
      if (all(!cell_actflag)) {
        
        write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
      
        thetamedian <- apply(theta_l, 2, median)
        thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
        thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
        
        write.csv(iz_l, file.path(sim_dir, paste0("iz_stage", l, ".csv")),row.names = FALSE)
        iz_cocluster_prob<-ArmListCoCluster(iz_l, Nsubgrp = N_subgrp, Narm=N_arm)
        
        final_efficacy_vec <- as.vector(t(cell_final_efficacy))
        write.csv(cell_final_efficacy, file.path(sim_dir, "final_efficacy.csv"),row.names = FALSE)
        
        final_futility_vec <- as.vector(t(cell_final_futility))
        write.csv(cell_final_futility, file.path(sim_dir, "final_futility.csv"),row.names = FALSE)
        write.csv(Status_Log, file.path(sim_dir, "status_log.csv"),row.names = FALSE)
        
        thetamedian_mat<-matrix(thetamedian,nrow = N_arm,ncol = N_subgrp,byrow=TRUE)
        
        cell_sample_size=as.matrix(table(Dat$i, Dat$k))
        
        break  
      }
      
      
    } else {
      
      
      # ---- Final analysis ----
      
      
      write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
      thetamedian <- apply(theta_l, 2, median)
      thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
      thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
 
      write.csv(iz_l, file.path(sim_dir, paste0("iz_stage", l, ".csv")),row.names = FALSE)
      iz_cocluster_prob<-ArmListCoCluster(iz_l, Nsubgrp = N_subgrp, Narm=N_arm)
      
      # draw final conclusion for active cells
      decision_l <- final_decision(
        theta_posterior = theta_l,
        cell_active = cell_actflag,
        cell_final_efficacy = cell_final_efficacy,
        cell_final_futility =cell_final_futility,
        futility_bound = bound_fut,
        efficacy_bound = bound_eff,
        futility_threshold = P_fut[l],
        efficacy_threshold = P_eff[l]
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
  return(
    list(
      sample_size=TotalSampleSize[l],
      cell_sample_size=as.vector(t(cell_sample_size)),
      theta=theta_l,
      thetamedian=thetamedian,
      thetaq750=thetaq750,
      thetaq250=thetaq250,
      iz_cocluster_prob=iz_cocluster_prob,
      final_efficacy=final_efficacy_vec,
      final_futility=final_futility_vec
    )
  )
}
