# -------------------------------------------------------------------
# Script: BHM_sim.R
# Purpose: Run simulation studies under multiple scenarios to evaluate
#          the performance of the BHM method on the same datasets
#          generated in BHARP_sim.R.
# Author: Xianglin Zhao
# Contents:
#   - Set root directory and load required packages and files
#   - Define hyperparameters
#   - For each scenario:
#       * Load datasets generated in BHARP_sim.R
#       * Run BHM in parallel
#       * Save posterior summaries
#       * Summarize results to log file
# Usage: Run BHARP_sim.R first to generate simulation datasets
#
# Dependencies:
#   R packages: here, rstan, posterior, parallel
#   Source file: BHM_model.R
#   External objects: project_root, scenario_list (from BHARP_sim.R)
# -------------------------------------------------------------------







# ---- Set root directory and load packages ----
setwd(project_root)     # Set project root manually
# Load packages after setting working directory
library(here);library(rstan);library(posterior);library(parallel)
#load source file
set_here()             # Set 'here' root path for relative file access
source(here("functions","BHM_model.R")) 


#----Setup----
N_sim      <- 500      # number of simulated trials
N_subgrp   <-10        # number of subgroups
N_arm      <-1         #number of arms

#----Hyperparameters----
BHM_hyperparam<-list(a_cell=5,b_cell=6,
                     mean_c0=0,pre_c0=2,
                     a_tau0= 4,b_tau0=4)


#----Simulation----
#scenario_list has been created in BHARP_sim.R
for (scenario in scenario_list) {
  scenario_dir <- here("results","simulation", scenario$name)
  #----1.Directory----
  # set working directory to S1...S9
  setwd(scenario_dir);
  # create a directory to save results 
  BHM_results_dir <- file.path(getwd(), "Output_BHM");if (!dir.exists(BHM_results_dir)) {dir.create(BHM_results_dir, showWarnings = FALSE, recursive = TRUE)}
  # data has been created by BHARP_sim.R 
  data_dir <- file.path(getwd(), "Data"); 
  # create log file
  log_file <- file.path(getwd(), paste0(scenario$name,"_BHM_sim_results.txt"))
  log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
  #----2.Extract settings----
  # extract settings 
  Theta_ <- scenario$Theta; ThetaVec<-as.vector(t(Theta_))
  cell_sizes <- scenario$cell_sizes
  
  
  # ----3.Conduct BHM----
  start<-Sys.time()
  cl <- makeCluster( detectCores() - 2,type="FORK")
  
  BHM_results <- parLapply(cl, 1:N_sim, function(v) {
    tryCatch({
      Dat<-read.csv(file.path(data_dir,paste0("DAT",v,".csv")))
      
      # get theta 
      BHM_result_v <- BHM_analysis(Data=Dat,
                                   n_subgrp=N_subgrp,
                                   n_arm=N_arm,
                                   hyperparam=BHM_hyperparam)
      
      thetamedian <- apply(BHM_result_v, 2, median)
      thetaq750   <- apply(BHM_result_v, 2, quantile, prob=0.750)
      thetaq250   <- apply(BHM_result_v, 2, quantile, prob=0.250)
      
      return(list(DAT_number=v,
                  theta=BHM_result_v, 
                  thetamedian=thetamedian, 
                  thetaq250=thetaq250,
                  thetaq750=thetaq750))
    }, error=function(e) {
      # if throw an error it will return a list with error message
      return(list(DAT_number=v,error=TRUE, message=e$message))
    })
  })
  
  stopCluster(cl);end<-Sys.time();round( end- start,2)
  
  BHM_good_results <- BHM_results[ sapply(BHM_results, function(x) is.null(x$error)) ]
  saveRDS(BHM_good_results, file = file.path(BHM_results_dir,"BHM_good_results_list.rds"))
  
  #collect results
  BHM_thetamedian <- do.call(rbind, lapply(BHM_good_results, "[[", "thetamedian"))
  
  
  log_message("====== BHM Results ======")
  
  # ----Root mean squared error----
  RMSE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      sqrt(mean(  (BHM_thetamedian[,x] - ThetaVec[x])^2 )  )
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("RMSE:")
    log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))
  }
  
  # ----Mean Absolute  Error----
  MAE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      mean( abs( (BHM_thetamedian[,x] - ThetaVec[x])  ) )  
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("Mean Abs Error:")
    log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))
  }
  
  
  # ----Variance of Estimates----
  VAR <- matrix(apply(BHM_thetamedian,2,var),nrow = N_arm,byrow=TRUE)
  {
    log_message("VAR of posterior median:")
    log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))
  }
  
}



