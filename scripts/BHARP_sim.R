# -------------------------------------------------------------------------
# BHARP Simulation Study
#
# This script reproduces the simulation study used to evaluate the
# Bayesian Hierarchical Adjustable Random Partition (BHARP) method.
#
# For each simulation scenario, the script:
#   1. generates independent normally distributed outcomes;
#   2. fits the BHARP model to each simulated dataset;
#   3. summarizes posterior estimates and subgroup partitions; 
#   4. reports operating characteristics.
#
# Before running:
#   - Set project_root to the local path of the cloned repository.
#   - Install the R packages listed under "Dependencies" below.
#
# Dependencies:
#   R packages:
#     Rcpp, RcppArmadillo, parallel, dplyr, purrr
#
#   Repository files:
#     functions/BHARP_theta.cpp
#     functions/BHARP_model.R
#     functions/helper_cocluster.R
#     scripts/simulation_config.R
#
# Output:
#   Generated datasets, posterior samples, summary statistics, and runtime
#   information are written to results/simulation/<scenario_name>/.
#
# Note:
#   The simulation uses FORK-based parallel processing and is intended for
#   macOS and Linux.
#
# -------------------------------------------------------------------------



# ---- Set root directory, load source files ----

project_root <- "YOURPATH/BHARP_release" #please modify this
if (!dir.exists(project_root)) {
  stop("ERROR: project_root does not exist. Please modify the project_root path at the top of the script.")
}


source(file.path(project_root, "functions", "helper_cocluster.R"))
Rcpp::sourceCpp(file.path(project_root, "functions", "BHARP_theta.cpp"))
source(file.path(project_root, "functions", "BHARP_model.R"))
source(file.path(project_root, "scripts", "simulation_config.R"))

# ---- Summarize BHARP simulation results ----
summarize_BHARP_results <- function( ThetaVec, N_arm, N_subgrp, good_results){
  All_thetamedian <- do.call(rbind, lapply(good_results, "[[", "thetamedian"))
  All_coCluster_list <- lapply(good_results, function(x) x$iz_cocluster_prob)
  All_thetaq250 <- do.call(rbind, lapply(good_results, "[[", "thetaq250"))
  All_thetaq750 <- do.call(rbind, lapply(good_results, "[[", "thetaq750"))
  All_theta_var <- do.call(rbind, lapply(good_results, "[[", "theta_var"))
  all_q_mode <- sapply(good_results, "[[", "q_mode")
  
  log_message("\n ====== Results ======")
  
  # Root mean squared error
  RMSE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      sqrt(mean(  (All_thetamedian[,x] - ThetaVec[x])^2 )  )
  ),nrow = N_arm,byrow=TRUE)
  log_message("RMSE:")
  log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))
  
  # Mean Absolute  Error
  MAE <- matrix(sapply(
    seq_len(N_arm*N_subgrp), function(x) 
      mean( abs( (All_thetamedian[,x] - ThetaVec[x])  ) )  
  ),nrow = N_arm,byrow=TRUE)
  log_message("Mean Abs Error:")
  log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))
  
  
  # ----Variance of Estimates----
  VAR <- matrix(apply(All_thetamedian,2,var),nrow = N_arm,byrow=TRUE)
  log_message("VAR of posterior median:")
  log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))

  
  
  
  # IQR
  IQR_mat <- All_thetaq750 - All_thetaq250
  IQR_mean <- matrix(colMeans(IQR_mat), nrow = N_arm, byrow = TRUE)
  log_message("Average posterior IQR of theta:")
  log_message(paste(capture.output(print(round(IQR_mean, 3))), collapse = "\n"))
  
  # posterior variance 
  PosteriorVar_mean <- matrix(colMeans(All_theta_var), nrow = N_arm, byrow = TRUE)
  log_message("Average posterior variance of theta (within simulations):")
  log_message(paste(capture.output(print(round(PosteriorVar_mean, 3))), collapse = "\n"))
  
  # distribution of q mode
  log_message("Posterior mode of q across simulated trials (relative frequency):")
  log_message(paste(capture.output(print(round(prop.table(table(all_q_mode)), 3))), collapse = "\n"))
  
  # Average Cocluster Probability
  Avg_coCluster_list <- vector("list", N_arm)  
  n_trials <- length(All_coCluster_list)      
  for (i in seq_len(N_arm)) {
    mats_for_arm_i <- lapply(All_coCluster_list, function(xx) xx[[i]])
    sum_mat <- Reduce("+", mats_for_arm_i)
    avg_mat <- round(sum_mat / n_trials,3)
    
    Avg_coCluster_list[[i]] <- avg_mat
  } 
  log_message("Coclustering Probability:")
  log_message(paste(capture.output(print(Avg_coCluster_list)), collapse = "\n"))
  
  return()
}


#---- Run simulations ----

for (s in seq_along(scenario_list)) {
  scenario <- scenario_list[[s]]
  set.seed(666)
 
  #---- 1.Directory ----
  
  scenario_dir <- file.path(
    project_root, 
    "results", 
    "simulation", 
    scenario$name
  )
  if (!dir.exists(scenario_dir)) {
    dir.create(scenario_dir, showWarnings = FALSE, recursive = TRUE)
  }
  results_dir <- file.path(
    scenario_dir, 
    "Output_BHARP"
  )
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  }
  data_dir <- file.path(
    scenario_dir, 
    "Data"
  )
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  }
  log_file <- file.path(
    scenario_dir, 
    paste0(scenario$name, "_BHARP_sim_results.txt")
  )
  log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
  #---- 2.Generate datasets ----
  Theta_ <- scenario$Theta
  ThetaVec <- as.vector(t(Theta_))
  cell_sizes <- scenario$cell_sizes
  
  Dat <- expand.grid(
    i = seq_len(N_arm),
    k = seq_len(N_subgrp)
  ) %>%
    purrr::pmap_dfr(function(i, k) {
      n <- cell_sizes[i, k]
      
      tibble::tibble(
        i = i,
        k = k,
        l = seq_len(n)
      )
    })
  
  log_message("====== Simulation setting ======")
  log_message(paste0("N_sim = ", N_sim, ", N_subgrp = ", N_subgrp, ", N_arm = ", N_arm))
  log_message(paste0("sample size = { " , paste(cell_sizes, collapse = ", ") , " }"))
  log_message(paste0("Theta = { " , paste(ThetaVec, collapse = ", ") , " }"))
  log_message("Hyperparameters:")
  log_message(paste(capture.output(str(hyperparam_list)), collapse = "\n"))
  
  
  for(v in 1:N_sim){
    Dat$Y <- rnorm(
      nrow(Dat),
      mean = Theta_[cbind(Dat$i, Dat$k)],
      sd = 1
    )
    write.csv(Dat, file.path(data_dir,paste0("DAT",v,".csv")),row.names=FALSE)
  }
  for (v in 1:N_sim) {
    f <- file.path(data_dir, paste0("DAT",v,".csv"))
    if (!file.exists(f)) cat("Missing file: ", f, "\n")
    if (file.info(f)$size == 0) cat("Empty file: ", f, "\n")
  }



  # ---- 3.Fit the BHARP model ----
  
  start<-Sys.time()
  
  cl <- parallel::makeCluster( max(1L, parallel::detectCores() - 2L), type="FORK")
  
  results <- parallel::parLapply(
    cl, 
    seq_len(N_sim), 
    function(v) {
      set.seed(66600000 + 1000 * s + v)
      tryCatch(
        {
          Dat<-read.csv(file.path(data_dir,paste0("DAT",v,".csv")))
          analysis_result_v <- Bharp_analysis(
            l = v,
            dataset=Dat,
            n_subgrp=N_subgrp,
            n_arm=N_arm,
            hyperparam=hyperparam_list, 
            results_dir=results_dir    
          )
          theta_v <- analysis_result_v$theta_posterior
          iz_v    <- analysis_result_v$iz_posterior
          iq_v    <- analysis_result_v$iq_posterior
      
          thetamedian <- apply(theta_v, 2, median)
          thetaq750   <- apply(theta_v, 2, quantile, prob=0.75)
          thetaq250   <- apply(theta_v, 2, quantile, prob=0.25)
          theta_var   <- apply(theta_v, 2, var)
      
          iz_cocluster_prob<-ArmListCoCluster(iz_v, Nsubgrp = N_subgrp, Narm=N_arm)
          q_tab   <- table(iq_v)
          q_mode  <- as.integer(names(q_tab)[which.max(q_tab)])
      
          list(
            DAT_number=v,
            theta=theta_v, 
            iz=iz_v,
            thetamedian=thetamedian, 
            thetaq250=thetaq250,
            thetaq750=thetaq750,
            theta_var=theta_var,
            iz_cocluster_prob=iz_cocluster_prob,
            q_mode=q_mode
          )
        }, 
        error=function(e) {
          list(DAT_number=v, error=TRUE, message=e$message)
        }
      )
    }
  )
  parallel::stopCluster(cl)

  scenario_time <- as.numeric(difftime(Sys.time(), start, units = "mins"))
  scenario_time <- round(scenario_time, 2)
  log_message(paste0("Total runtime: ", scenario_time, " minutes"))

  successful <- vapply(
    results,
    function(result) is.null(result$error),
    logical(1)
  )
  
  BHARP_good_results <- results[successful]
  BHARP_bad_results <- results[!successful]
  
  log_message( paste0( "Successful simulations: ", length(BHARP_good_results), " / ", N_sim ))
  
  log_message( paste0( "Failed simulations: ", length(BHARP_bad_results), " / ",N_sim ))

  if (length(BHARP_good_results) == 0L) {
    stop("All simulations failed. Check Bharp_analysis and C++ source.")
  }

  saveRDS(
    BHARP_good_results, 
    file = file.path(results_dir,"BHARP_good_results_list.rds")
  )


# ---- 4.Summarize operating characteristics ----
  summarize_BHARP_results( 
    ThetaVec, 
    N_arm, 
    N_subgrp, 
    BHARP_good_results
  )
  

}
