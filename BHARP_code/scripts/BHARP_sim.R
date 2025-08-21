# -------------------------------------------------------------------
# Script: BHARP_sim.R
# Purpose: Run simulation studies under multiple scenarios to evaluate
#          the performance of the BHARP method.
# Author: Xianglin Zhao
# Contents:
#   - Set root directory and load required packages and C++ files
#   - Define hyperparameters and simulation scenarios
#   - For each scenario:
#       * Generate data
#       * Run BHARP in parallel
#       * Save posterior summaries
#       * Summarize results to log file
# Usage: Modify `project_root` to your local directory
#
# Dependencies:
#   R packages: here, Rcpp, parallel, tidyr, dplyr, purrr
#   C++ file: BHARP.cpp
#   Helper functions: BHARP_model.R helper_metric.R
# -------------------------------------------------------------------


# ---- Set root directory and load packages ----
# Modify the path below to your actual local path
project_root <- "/YOURNAME/YOURPROJECT/BHARP_code" #please modify this
setwd(project_root)
# Load packages after setting working directory
library(here);library(Rcpp);library(parallel);library(tidyr);library(dplyr);library(purrr)

#load source file
set_here()
sourceCpp(here("functions","BHARP.cpp")) 
source(here("functions","BHARP_model.R")) 
source(here("functions","helper_metric.R")) 

#----Setup----
N_sim      <- 500      # number of simulated trials
N_subgrp   <-10        # number of subgroups
N_arm      <-1         #number of arms

Dat <- expand.grid(i = 1:N_arm, k = 1:N_subgrp) %>%
  pmap_dfr(function(i, k) {
    n <- cell_sizes[i, k]
    tibble(i = i, k = k, l = 1:n)
  })

#----Hyperparameters----
hyperparam_list<-list(c_beta=rep(0.0,N_arm), p_beta=rep(2.0,N_arm),
                      a_cell=5,b_cell=6,
                      a_between=4,b_between=0.8*5,
                      a_within=70,b_within=0.10^2*71)
#----Scenario----
scenario_list <- list(
  S1 = list(
    name = "S1",
    Theta = matrix(c(rep(0.00,10)), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S2 = list(
    name = "S2",
    Theta = matrix(c(rep(0.00,7), rep(1.3,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes= matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S3 = list(
    name = "S3",
    Theta = matrix(c(rep(0.00,7), rep(0.65,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S4 = list(
    name = "S4",
    Theta = matrix(c(rep(0.00,7), rep(0.65,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(c(rep(35,7),rep(70,3)), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S5 = list(
    name = "S5",
    Theta = matrix(c(rep(0.00,5), rep(0.65,5) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S6 = list(
    name = "S6",
    Theta = matrix(c(rep(0.00,4), rep(0.65,3), rep(1.30,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S7 = list(
    name = "S7",
    Theta = matrix(c(rep(0.00,4), rep(0.65,3), rep(1.30,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(c(rep(42,4),rep(84,3),rep(56,3)), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S8 = list(
    name = "S8",
    Theta = matrix(c(rep(0.00,7), rep(0.65,2), rep(1.30,1) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S9 = list(
    name = "S9",
    Theta = matrix(c(rep(0.00,7), rep(0.65,2), rep(1.30,1) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(c(rep(35,7),rep(70,2),84), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) )
)


#----Simulation----

for (scenario in scenario_list) {
  set.seed("666")
  scenario_dir <- here("results","simulation", scenario$name)
  if (!dir.exists(scenario_dir)) {dir.create(scenario_dir, showWarnings = FALSE, recursive = TRUE)}
  #----1.Directory----
  # set working directory to S1...S9
  setwd(scenario_dir);
  # create a directory to save results 
  results_dir <- file.path(getwd(), "Output_BHARP"); if (!dir.exists(results_dir)) {dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)}
  # create directory to save data
  data_dir <- file.path(getwd(), "Data"); if (!dir.exists(data_dir)) {dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)}
  # create log file
  log_file <- file.path(getwd(), paste0(scenario$name,"_BHARP_sim_results.txt"))
  log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
  #----2.Data generating----
  # extract settings 
 
  Theta_ <- scenario$Theta; ThetaVec<-as.vector(t(Theta_))
  cell_sizes <- scenario$cell_sizes
  
  log_message("====== Simulation setting ======")
  log_message(paste0("N_sim = ", N_sim, ", N_subgrp = ", N_subgrp, ", N_arm = ", N_arm))
  log_message(paste0("sample size = { " , paste(cell_sizes, collapse = ", ") , " }"))
  log_message(paste0("Theta = { " , paste(ThetaVec, collapse = ", ") , " }"))
  log_message("Hyperparameters:")
  log_message(paste(capture.output(str(hyperparam_list)), collapse = "\n"))
  
  # generate and save datasets  
  for(v in 1:N_sim){
    Dat$Y<-rnorm(nrow(Dat), mean=Theta_[cbind(Dat$i, Dat$k)], sd=1)
    write.csv(Dat, file.path(data_dir,paste0("DAT",v,".csv")),row.names=FALSE)
  }
  # check availability
  for (v in 1:N_sim) {
    f <- file.path(data_dir, paste0("DAT",v,".csv"))
    if (!file.exists(f)) cat("Missing file: ", f, "\n")
    if (file.info(f)$size == 0) cat("Empty file: ", f, "\n")
  }



  # ----3.Conduct BHARP----
  set.seed(1);start<-Sys.time()
  cl <- makeCluster( detectCores() - 2,type="FORK")
  results <- parLapply(cl, 1:N_sim, function(v) {
    tryCatch({
      Dat<-read.csv(file.path(data_dir,paste0("DAT",v,".csv")))
      # get theta and iz
      analysis_result_v <- Bharp_analysis(v,
                                          dataset=Dat,
                                          n_subgrp=N_subgrp,
                                          n_arm=N_arm,
                                          hyperparam=hyperparam_list, 
                                          results_dir    )
      theta_v <- analysis_result_v$theta_posterior
      iz_v    <- analysis_result_v$iz_posterior
      
      thetamedian <- apply(theta_v, 2, median)
      thetaq750   <- apply(theta_v, 2, quantile, prob=0.75)
      thetaq250   <- apply(theta_v, 2, quantile, prob=0.25)
      iz_cocluster_prob<-ArmListCoCluster(iz_v, Nsubgrp = N_subgrp, Narm=N_arm)
      return(list(DAT_number=v,
                  theta=theta_v, iz=iz_v,
                  thetamedian=thetamedian, 
                  thetaq250=thetaq250,
                  thetaq750=thetaq750,
                  iz_cocluster_prob=iz_cocluster_prob))
    }, error=function(e) {
      # if throw an error it will return a list with error message
      return(list(DAT_number=v,error=TRUE, message=e$message))
    })
  })
  stopCluster(cl);round(Sys.time() - start,2)


# results is a list of N_sims 
# for cycles with error, results[[v]]$error=TRUE, otherwise results[[v]]$thetamedian etc have contents
good_results <- results[ sapply(results, function(x) is.null(x$error)) ]
saveRDS(good_results, file = file.path(results_dir,"BHARP_good_results_list.rds"))


# ----4.Calculate OCs----
summarize_BHARP_results( ThetaVec, N_arm, N_subgrp, N_sim, BHARP_good_results)
  

}
