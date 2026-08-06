# -------------------------------------------------------------------------
# BHARP Adaptive Trial Design Evaluation
#
# This script evaluates the operating characteristics of the BHARP
# adaptive trial design through repeated trial simulations.
#
# For each simulated trial, the script:
#   1. recruits patients according to the prespecified sample sizes;
#   2. performs interim and final BHARP analyses;
#   3. applies the futility and efficacy decision rules; 
#   4. records estimation and design operating characteristics.
#
# Before running:
#   - Set project_root to the local path of the cloned repository.
#   - Install the required R packages.
#
# Dependencies:
#   R packages:
#     Rcpp, RcppArmadillo
#
#   The parallel package is included with R.
#
#   Repository files:
#     functions/BHARP_theta.cpp
#     functions/BHARP_model.R
#     functions/helper_cocluster.R
#     functions/helper_trial.R
#     scripts/trial_design_config.R
#
# Output:
#   Trial-level data and posterior results are written to:
#     results/trial_design/trial_BHARP/
#   Summary operating characteristics are written to:
#     results/trial_design/BHARP_trial_log.txt
#
# Note:
#   The simulation uses FORK-based parallel processing and is intended for
#   macOS and Linux.
# -------------------------------------------------------------------------


# ---- Set root directory, load source files ----

project_root <- "YOURPATH/BHARP_release" #please modify this
if (!dir.exists(project_root)) {
  stop("ERROR: project_root does not exist. Please modify the project_root path at the top of the script.")
}

results_dir  <- file.path(
  project_root, 
  "results",
  "trial_design",
  "trial_BHARP"
)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
}


sourceCpp(file.path(project_root, "functions", "BHARP_theta.cpp"))
source(file.path(project_root, "functions", "helper_cocluster.R"))
source(file.path(project_root, "functions", "helper_trial.R"))
source(file.path(project_root, "functions", "BHARP_model.R"))
source(file.path(project_root, "scripts", "trial_design_config.R"))








# ---- Creating simulation log ----
log_file <- file.path(dirname(results_dir), "BHARP_trial_log.txt")
log_message <- function(msg) { 
  cat(paste0( msg, "\n"), file = log_file, append = TRUE)
}
log_message("====== Simulation setting ======")
log_message(paste0("N_sim = ", N_sim, ", N_subgrp = ", N_subgrp,
                     ", N_arm = ", N_arm, ", L = ", L))
log_message(paste0("x_fut = ", x_fut, ", x_eff = ", x_eff))
log_message(paste0("P_fut = ", paste(P_fut, collapse = ", "),
                     ", P_eff = ", paste(P_eff, collapse = ", ")))
log_message(paste0("TotalSampleSize = { " , paste(TotalSampleSize, collapse = ", ") , " }"))
log_message("Theta matrix (true cell means):")
log_message(paste(capture.output(print(Theta)), collapse = "\n"))
log_message("Hyperparameters:")
log_message(paste(capture.output(str(hyperparam_list)), collapse = "\n"))
  

  




# ---- Run simulated trials ----

start <- Sys.time() 

cl <- parallel::makeCluster(max(1L,parallel::detectCores() - 2L), type="FORK")

results <- parallel::parLapply(
  cl, 
  seq_len(N_sim), 
  function(v) {
    trial_seed <- 6660000 + v
    set.seed(trial_seed)
    tryCatch(
      {
        out <- simulate_one_trial(
          trial_index = v, 
          N_subgrp = N_subgrp, 
          N_arm = N_arm,
          L=L,
          TotalSampleSize = TotalSampleSize,
          Theta=Theta,
          hyperparam_list = hyperparam_list,
          bound_fut = x_fut, 
          bound_eff = x_eff, 
          P_fut=P_fut, 
          P_eff=P_eff,
          base_dir=results_dir
        )
        out$trial_index <- v
        out$seed <- trial_seed
        out
      }, 
      error=function(e) {
        list(error=TRUE,trial_index = v, message=e$message)
      }
    )
  }
)

parallel::stopCluster(cl); 

scenario_time <- as.numeric(difftime(Sys.time(), start, units = "mins"))
scenario_time <- round(scenario_time, 2)
log_message(paste0("Total runtime: ", scenario_time, " minutes"))





successful <- vapply(
  results,
  function(result) is.null(result$error),
  logical(1)
)

good_results <- results[successful]
bad_results <- results[!successful]

n_good <- length(good_results)
if (n_good == 0L) {
  stop( "All simulated trials failed. Check the trial-level error messages." )
}


saveRDS(good_results, file =file.path(results_dir,"BHARP_good_results_list.rds"))
log_message(paste0("Failed trials: ", length(bad_results), " / ", N_sim))


All_thetamedian <- do.call(rbind, lapply(good_results, "[[", "thetamedian"))
All_thetaq750 <- do.call(rbind,lapply(good_results, "[[", "thetaq750"))
All_thetaq250 <- do.call(rbind,lapply(good_results, "[[", "thetaq250"))

All_thetapostvar <- do.call(rbind, lapply(good_results, function(x) {
  apply(x$theta, 2, var)
}))

All_efficacy <- do.call(rbind,lapply(good_results, "[[", "final_efficacy"))
All_cell_sample_size <- do.call(rbind,lapply(good_results, "[[", "cell_sample_size"))
All_coCluster_list <- lapply(good_results, function(x) x$iz_cocluster_prob)




# ---- Calculate operating characteristics ----

#Mean Absolute Error
MAE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    mean( abs( (All_thetamedian[,x] - ThetaVec[x])  ) )  
),nrow = N_arm,ncol=N_subgrp,byrow = TRUE)

log_message("Mean Abs Error:")
log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))




# Variance of posterior median estimates across simulated trials
VAR<- matrix(apply(All_thetamedian,2,var),nrow = N_arm,ncol=N_subgrp,byrow = TRUE)
log_message("VAR of posterior median:")
log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))


# Root mean squared error
RMSE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    sqrt(mean(  (All_thetamedian[,x] - ThetaVec[x])^2 )  )
),nrow = N_arm,ncol=N_subgrp,byrow = TRUE)

log_message("RMSE:")
log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))







# Average posterior Inter-Quantile Ranges
IQR<-matrix(colMeans(All_thetaq750-All_thetaq250),nrow = N_arm,ncol=N_subgrp,byrow = TRUE)

log_message("Inter-Quantile Ranges:")
log_message(paste(capture.output(print(round(IQR,3))), collapse = "\n"))


# posterior variance
PostVar<-matrix(colMeans(All_thetapostvar),nrow = N_arm,ncol=N_subgrp,byrow = TRUE)
log_message("Posterior Variance:")
log_message(paste(capture.output(print(round(PostVar,3))), collapse = "\n"))








#Expected Cell Sample Size
CellSampleSize<-matrix(colMeans(All_cell_sample_size),nrow = N_arm,ncol=N_subgrp,byrow = TRUE)
log_message("Expected Cell Sample Size:")
log_message(paste(capture.output(print(round(CellSampleSize))), collapse = "\n"))



#Global False Positive Rate
true_effective <- ThetaVec>x_eff
true_effective_mat <- matrix(rep(true_effective,each=n_good),nrow = n_good,ncol = N_arm*N_subgrp)
fp_mat=(!true_effective_mat)& All_efficacy
GFPR<-mean(apply(fp_mat, 1, any))
log_message("Global False Positive Rate:")
log_message(GFPR)





#Generalized Power by subgrp
GP_grp<-mean(sapply(seq_len(n_good), function(v) 
trial_success_subgroup(All_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by subgrp:")
log_message(GP_grp)

# Generalized Power by arm
GP_arm<-mean(sapply(seq_len(n_good), function(v) 
trial_success_arm(All_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by arm:")
log_message(GP_arm)


# Average Cocluster Probability

Avg_coCluster_list <- vector("list", N_arm)  # output list
n_trials <- length(All_coCluster_list)      # valid trials

for (i in seq_len(N_arm)) {
  mats_for_arm_i <- lapply(All_coCluster_list, function(xx) xx[[i]])
  sum_mat <- Reduce("+", mats_for_arm_i)
  avg_mat <- round(sum_mat / n_trials,3)
  Avg_coCluster_list[[i]] <- avg_mat
} 

log_message("Coclustering Probability:")
log_message(paste(capture.output(print(Avg_coCluster_list)), collapse = "\n"))






