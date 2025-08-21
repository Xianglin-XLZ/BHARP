# File: helper_metric.R
# Purpose: Define utility functions to evaluate generalized power after simulation
# Author:Xianglin Zhao
# Contents:
# - trial_success_subgroup(): Assess power by detecting all effective subgroups
# - trial_success_arm(): Assess power by detecting all effective arms
# - summarize_BHARP_results(): log RMSE MAE VAR and cocluster probability of BHARP  

# Dependency: Requires `log_message()` and `integrated_RMSE()` to be defined in the global environment; input `good_results` should be a list of valid simulation outputs
# These functions are used *after* all simulations are completed to evaluate estimation accuracy
# and overall success from different perspectives (e.g., by subgroup or intervention).







# Generalized power by subgroup: declare success only if all truly sensitive subgroups are correctly detected
# For each truly effective subgroup, check if at least one effective arm is identified

trial_success_subgroup<-function(efficacy_conclusion, true_effective, N_subgrp, N_arm){
  efficacy_conclusion_Mat <- matrix(efficacy_conclusion, ncol=N_subgrp,nrow=N_arm,byrow=TRUE)
  true_effective_Mat <- matrix(true_effective,           ncol=N_subgrp,nrow=N_arm,byrow=TRUE)
  
  true_sensitive_subgroups <- which(colSums(true_effective_Mat) > 0)
  successful_detect <- rep(NA,N_subgrp)
  for (k in true_sensitive_subgroups) {
    #correctly detect any useful arm
    successful_detect[k]  <- any(true_effective_Mat[,k] & efficacy_conclusion_Mat[,k])
  }
  return(all(successful_detect,na.rm=TRUE))
}


# Generalized power by arm: declare success only if all truly effective interventions are correctly detected
# For each truly effective arm, check if at least one effective subgroup is identified
trial_success_arm<-function(efficacy_conclusion, true_effective, N_subgrp, N_arm){
  efficacy_conclusion_Mat <- matrix(efficacy_conclusion, ncol=N_subgrp,nrow=N_arm,byrow=TRUE)
  true_effective_Mat <- matrix(true_effective,           ncol=N_subgrp,nrow=N_arm,byrow=TRUE)
  
  true_useful_arms <- which(rowSums(true_effective_Mat) > 0)
  successful_detect <- rep(NA,N_arm)
  for (i in true_useful_arms) {
    #correctly detect any sensitive subgroup
    successful_detect[i]  <- any(true_effective_Mat[i,] & efficacy_conclusion_Mat[i,])
  }
  return(all(successful_detect,na.rm=TRUE))
}


  

# log RMSE MAE VAR and cocluster probability of BHARP in simulation study 
summarize_BHARP_results <- function( ThetaVec, N_arm, N_subgrp, N_sim, good_results){
  
  # collect results 
  All_thetamedian <- do.call(rbind, lapply(good_results, "[[", "thetamedian"))
  All_coCluster_list <- lapply(good_results, function(x) x$iz_cocluster_prob)

  
  log_message("\n ====== Results ======")
  
  # ----Root mean squared error----
  RMSE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      sqrt(mean(  (All_thetamedian[,x] - ThetaVec[x])^2 )  )
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("RMSE:")
    log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))
  }

  # ----Mean Absolute  Error----
  MAE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      mean( abs( (All_thetamedian[,x] - ThetaVec[x])  ) )  
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("Mean Abs Error:")
    log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))
  }
  
  
  # ----Variance of Estimates----
  VAR <- matrix(apply(All_thetamedian,2,var),nrow = N_arm,byrow=TRUE)
  {
    log_message("VAR of posterior median:")
    log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))
  }
  

  # ----Average Cocluster Probability----
  Avg_coCluster_list <- vector("list", N_arm)  # output list
  n_trials <- length(All_coCluster_list)      # valid trials
  
  #  i=1..Narm
  for (i in seq_len(N_arm)) {
    mats_for_arm_i <- lapply(All_coCluster_list, function(xx) xx[[i]])
    sum_mat <- Reduce("+", mats_for_arm_i)
    avg_mat <- round(sum_mat / n_trials,3)
    
    Avg_coCluster_list[[i]] <- avg_mat
  } 
  {
    log_message("Coclustering Probability:")
    log_message(paste(capture.output(print(Avg_coCluster_list)), collapse = "\n"))
  }
  return()
}
