# -------------------------------------------------------------------
# Script: BLAST_sim.R
# Purpose: Run simulation studies under multiple scenarios to evaluate
#          the performance of the BLAST method on the same datasets
#          generated in BHARP_sim.R.
# Author: Xianglin Zhao
# Contents:
#   - Set root directory and load required packages and files
#   - Define hyperparameters
#   - For each scenario:
#       * Load datasets generated in BHARP_sim.R
#       * Run BLAST in parallel
#       * Save posterior summaries
#       * Summarize results to log file
# Usage: Run BHARP_sim.R first to generate simulation datasets
#
# Dependencies:
#   R packages: here, rstan, posterior, parallel
#   Source file: BLAST_model.R
#   External objects: project_root, scenario_list (from BHARP_sim.R)
# -------------------------------------------------------------------



# ---- Set root directory and load packages ----
setwd(project_root)
# Load packages after setting working directory
library(here);library(rstan);library(posterior);library(parallel)
#load source file
set_here()
source(here("functions","BLAST_model.R")) 

#----Setup----
N_sim      <- 500      # number of simulated trials
N_subgrp   <-10        # number of subgroups
N_arm      <-1         #number of arms

#----Hyperparameters----
BLAST_hyperparam<-list(c_beta=0,p_beta=2,
                       a_cell=5,b_cell=6,
                       a_between= 4, b_between=4,
                       a_within=70, b_within=0.71)



#----Simulation----
#scenario_list has been created in BHARP_sim.R
for (scenario in scenario_list) {
  scenario_dir <- here("results","simulation", scenario$name)
  #----1.Directory----
  # set working directory to S1...S9
  setwd(scenario_dir);
  # create a directory to save results 
  BLAST_results_dir <- file.path(getwd(), "Output_BLAST");if (!dir.exists(BLAST_results_dir)) {dir.create(BLAST_results_dir, showWarnings = FALSE, recursive = TRUE)}
  # data has been created by BHARP_sim.R 
  data_dir <- file.path(getwd(), "Data"); 
  # create log file
  log_file <- file.path(getwd(), paste0(scenario$name,"_BLAST_sim_results.txt"))
  log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
  #----2.Extract settings----
  # extract settings 
  Theta_ <- scenario$Theta; ThetaVec<-as.vector(t(Theta_))
  cell_sizes <- scenario$cell_sizes
  





  # ----3. Conduct BLAST----
  start<-Sys.time()
  cl <- makeCluster( detectCores() - 2,type="FORK")
  
  BLAST_results <- parLapply(cl, 1:N_sim, function(v) {
    tryCatch({
      Dat<-read.csv(file.path(data_dir,paste0("DAT",v,".csv")))
      
      BLAST<-BLAST_analysis(Data=Dat,
                            n_subgrp=N_subgrp,
                            hyperparam=BLAST_hyperparam)
      
      # get theta 
      BLAST_result_v <- BLAST$theta
      BLAST_q_v<-BLAST$q
      
      
      thetamedian <- apply(BLAST_result_v, 2, median)
      thetaq750   <- apply(BLAST_result_v, 2, quantile, prob=0.750)
      thetaq250   <- apply(BLAST_result_v, 2, quantile, prob=0.250)
      
      
      
      return(list(DAT_number=v,
                  q=BLAST_q_v,
                  theta=BLAST_result_v, 
                  thetamedian=thetamedian, 
                  thetaq250=thetaq250,
                  thetaq750=thetaq750))
    }, error=function(e) {
      # if throw an error it will return a list with error message
      return(list(DAT_number=v,error=TRUE, message=e$message))
    })
  })
  
  stopCluster(cl);stop<-Sys.time();round(stop - start,2)
  
  BLAST_good_results <- BLAST_results[ sapply(BLAST_results, function(x) is.null(x$error)) ]
  saveRDS(BLAST_good_results, file = file.path(BLAST_results_dir,"BLAST_good_results_list.rds"))
  
  
  BLAST_thetamedian <- do.call(rbind, lapply(BLAST_good_results, "[[", "thetamedian"))
  BLAST_thetaq750 <- do.call(rbind,lapply(BLAST_good_results, "[[", "thetaq750"))
  BLAST_thetaq250 <- do.call(rbind,lapply(BLAST_good_results, "[[", "thetaq250"))
  BLAST_q<-do.call(c,lapply(BLAST_good_results, "[[", "q"))
  
  
  log_message("====== BLAST Results ======")
  # ----Root mean squared error----
  RMSE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      sqrt(mean(  (BLAST_thetamedian[,x] - ThetaVec[x])^2 )  )
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("RMSE:")
    log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))
  }
  
  # ----Mean Absolute  Error----
  MAE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      mean( abs( (BLAST_thetamedian[,x] - ThetaVec[x])  ) )  
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("Mean Abs Error:")
    log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))
  }
  
  
  # ----Variance of Estimates----
  VAR <- matrix(apply(BLAST_thetamedian,2,var),nrow = N_arm,byrow=TRUE)
  {
    log_message("VAR of posterior median:")
    log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))
  }

}