# -------------------------------------------------------------------
# Script: BART_sim.R
# Purpose: Run simulation studies under multiple scenarios to evaluate
#          the performance of the BART method on the same datasets
#          generated in BHARP_sim.R.
# Author: Xianglin Zhao
# Contents:
#   - Set root directory and load required packages and files
#   - Define hyperparameters
#   - For each scenario:
#       * Load datasets generated in BHARP_sim.R
#       * Run BART in parallel
#       * Save posterior summaries
#       * Summarize results to log file
# Usage: Run BHARP_sim.R first to generate simulation datasets
#
# Dependencies:
#   R packages: here, BART, posterior, parallel
#   Source file: BART_model.R
#   External objects: project_root, scenario_list (from BHARP_sim.R)
# -------------------------------------------------------------------


# ---- Set root directory and load packages ----
setwd(project_root)
# Load packages after setting working directory
library(here);library(BART);library(posterior);library(parallel)
#load source file
set_here()
source(here("functions","BART_model.R")) 


#----Setup----
N_sim      <- 500      # number of simulated trials
N_subgrp   <-10        # number of subgroups
N_arm      <-1         #number of arms






#----Simulation----
#scenario_list has been created in BHARP_sim.R
for (scenario in scenario_list) {
  scenario_dir <- here("results","simulation", scenario$name)
  #----1.Directory----
  # set working directory to S1...S9
  setwd(scenario_dir);
  # create a directory to save results 
  BART_results_dir <- file.path(getwd(), "Output_BART");if (!dir.exists(BART_results_dir)) {dir.create(BART_results_dir, showWarnings = FALSE, recursive = TRUE)}
  # data has been created by BHARP_sim.R 
  data_dir <- file.path(getwd(), "Data"); 
  # create log file
  log_file <- file.path(getwd(), paste0(scenario$name,"_BART_sim_results.txt"))
  log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
  #----2.Extract settings----
  # extract settings 
  Theta_ <- scenario$Theta; ThetaVec<-as.vector(t(Theta_))
  cell_sizes <- scenario$cell_sizes
  
  
  # ----3. Conduct BART----
  start<-Sys.time()
  cl <- makeCluster( detectCores() - 2,type="FORK")
  
  BART_results <- parLapply(cl, 1:N_sim, function(v) {
    tryCatch({
      
      Dat<-read.csv(file.path(data_dir,paste0("DAT",v,".csv")))
      if(N_subgrp>1){
        Dat$k <- factor(Dat$k)
        k_dummies <- model.matrix(~ k - 1, data = Dat)
        Dat<-cbind(Dat,k_dummies)
      }
      
      one_BART_list<-list()
      
      for(i in 1:N_arm){
        Dat_i<-Dat[Dat$i==i,]
        X_i<-as.matrix(Dat_i[,paste0("k",1:N_subgrp)])
        Y_i<-as.matrix(Dat_i[,"Y"])
        
        bart_model <- wbart(x.train = X_i, 
                            y.train = Y_i,ndpost=2000)   
        yhat_matrix <- bart_model$yhat.train #  columns stand for individual
        
        for(k in 1:N_subgrp){
          idx_k <- which(Dat_i$k == k)
          yhat_k <- yhat_matrix[, idx_k, drop = FALSE]  # ndpost × n_k
          cell_mean_samples <- rowMeans(yhat_k)         # ndpost × 1
          one_BART_list[[paste0("Arm", i, "Grp", k)]] <- cell_mean_samples
        }
      }
      
      # get theta 
      BART_result_v <- do.call(cbind, one_BART_list)
      colnames(BART_result_v) <- names(one_BART_list)
      
      thetamedian <- apply(BART_result_v, 2, median)
      thetaq750   <- apply(BART_result_v, 2, quantile, prob=0.750)
      thetaq250   <- apply(BART_result_v, 2, quantile, prob=0.250)
      
      return(list(DAT_number=v,
                  theta=BART_result_v, 
                  thetamedian=thetamedian, 
                  thetaq250=thetaq250,
                  thetaq750=thetaq750))
    }, error=function(e) {
      # if throw an error it will return a list with error message
      return(list(DAT_number=v,error=TRUE, message=e$message))
    })
  })
  
  stopCluster(cl);end<-Sys.time();round(end-start,1)
  
  BART_good_results <- BART_results[ sapply(BART_results, function(x) is.null(x$error)) ]
  saveRDS(BART_good_results, file = file.path(BART_results_dir,"BART_good_results_list.rds"))
  
  
  BART_thetamedian <- do.call(rbind, lapply(BART_good_results, "[[", "thetamedian"))
  BART_thetaq750 <- do.call(rbind,lapply(BART_good_results, "[[", "thetaq750"))
  BART_thetaq250 <- do.call(rbind,lapply(BART_good_results, "[[", "thetaq250"))
  
  
  
  log_message("====== BART Results ======")
  # ----Root mean squared error----
  RMSE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      sqrt(mean(  (BART_thetamedian[,x] - ThetaVec[x])^2 )  )
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("RMSE:")
    log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))
  }
  
  # ----Mean Absolute  Error----
  MAE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      mean( abs( (BART_thetamedian[,x] - ThetaVec[x])  ) )  
  ),nrow = N_arm,byrow=TRUE)
  {
    log_message("Mean Abs Error:")
    log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))
  }
  
  
  # ----Variance of Estimates----
  VAR <- matrix(apply(BART_thetamedian,2,var),nrow = N_arm,byrow=TRUE)
  {
    log_message("VAR of posterior median:")
    log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))
  }
}





