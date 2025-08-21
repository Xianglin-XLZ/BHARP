# -------------------------------------------------------------------
# Script: IND_sim.R
# Purpose: Run simulation studies under multiple scenarios to evaluate
#          the performance of the IND method on the same datasets
#          generated in BHARP_sim.R.
# Author: Xianglin Zhao
# Contents:
#   - Set root directory and load required packages and files
#   - Define hyperparameters
#   - For each scenario:
#       * Load datasets generated in BHARP_sim.R
#       * Run IND in parallel
#       * Save posterior summaries
#       * Summarize results to log file
# Usage: Run BHARP_sim.R first to generate simulation datasets
#
# Dependencies:
#   R packages: here, rstan, posterior, parallel
#   Source file: IND_model.R
#   External objects: project_root, scenario_list (from BHARP_sim.R)
# -------------------------------------------------------------------

# ---- Set root directory and load packages ----
setwd(project_root)
# Load packages after setting working directory
library(here);library(rstan);library(posterior);library(parallel)
#load source file
set_here()
source(here("functions","IND_model.R")) 

#----Setup----
N_sim      <- 500      # number of simulated trials
N_subgrp   <-10        # number of subgroups
N_arm      <-1         #number of arms

#----Hyperparameters----
IND_hyperparam<-list(a_cell=5,b_cell=6,
                     c_thetaind=0,p_thetaind=2)



#----Simulation----
#scenario_list has been created in BHARP_sim.R
for (scenario in scenario_list) {
  scenario_dir <- here("results","simulation", scenario$name)
  #----1.Directory----
  # set working directory to S1...S9
  setwd(scenario_dir);
  # create a directory to save results 
  IND_results_dir <- file.path(getwd(), "Output_IND");if (!dir.exists(IND_results_dir)) {dir.create(IND_results_dir, showWarnings = FALSE, recursive = TRUE)}
  # data has been created by BHARP_sim.R 
  data_dir <- file.path(getwd(), "Data"); 
  # create log file
  log_file <- file.path(getwd(), paste0(scenario$name,"_IND_sim_results.txt"))
  log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
  #----2.Extract settings----
  # extract settings 
  Theta_ <- scenario$Theta; ThetaVec<-as.vector(t(Theta_))
  cell_sizes <- scenario$cell_sizes

  
  # ----3. Conduct IND----
  start<-Sys.time()
  cl <- makeCluster( detectCores() - 2,type="FORK")
  
  IND_results <- parLapply(cl, 1:N_sim, function(v) {
    tryCatch({
      Dat<-read.csv(file.path(data_dir,paste0("DAT",v,".csv")))
      
      # get theta 
      IND_result_v <- IND_analysis(Data=Dat,
                                   n_subgrp=N_subgrp,
                                   n_arm=N_arm,
                                   hyperparam=IND_hyperparam)
      
      thetamedian <- apply(IND_result_v, 2, median)
      thetaq750   <- apply(IND_result_v, 2, quantile, prob=0.750)
      thetaq250   <- apply(IND_result_v, 2, quantile, prob=0.250)
      return(list(DAT_number=v,
                  theta=IND_result_v, 
                  thetamedian=thetamedian, 
                  thetaq250=thetaq250,
                  thetaq750=thetaq750))
    }, error=function(e) {
      # if throw an error it will return a list with error message
      return(list(DAT_number=v,error=TRUE, message=e$message))
    })
  })
  
  stopCluster(cl);end<-Sys.time();round(end - start,2)
  
  IND_good_results <- IND_results[ sapply(IND_results, function(x) is.null(x$error)) ]
  saveRDS(IND_good_results, file = file.path(IND_results_dir,"IND_good_results_list.rds"))
  
  #collect results
  IND_thetamedian <- do.call(rbind, lapply(IND_good_results, "[[", "thetamedian"))

  
  
  log_message("====== IND Results ======")
  # ----Root mean squared error----
  RMSE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      sqrt(mean(  (IND_thetamedian[,x] - ThetaVec[x])^2 )  )
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("RMSE:")
    log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))
  }
  
  # ----Mean Absolute  Error----
  MAE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      mean( abs( (IND_thetamedian[,x] - ThetaVec[x])  ) )  
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("Mean Abs Error:")
    log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))
  }
  
  
  # ----Variance of Estimates----
  VAR <- matrix(apply(IND_thetamedian,2,var),nrow = N_arm,byrow=TRUE)
  {
    log_message("VAR of posterior median:")
    log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))
  }

}






