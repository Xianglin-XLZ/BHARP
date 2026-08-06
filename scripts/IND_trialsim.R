# -------------------------------------------------------------------------
# IND Adaptive Trial Design Evaluation
#
# This script evaluates the operating characteristics of the independent
# analysis (IND) adaptive trial design through repeated trial simulations.
#
# For each trial, the script performs stage-wise recruitment, interim and
# final IND analyses, and the shared futility and efficacy decisions.
#
# Before running:
#   - Set project_root to the local path of the cloned repository.
#   - Install the rstan package.
#   - Review the settings in scripts/trial_design_config.R.
#
# Dependencies:
#   R package:
#     rstan
#
#   Repository files:
#     functions/IND_model.R
#     functions/helper_trial.R
#     scripts/trial_design_config.R
#
# Output:
#   Trial-level results:
#     results/trial_design/trial_IND/
#
#   Summary operating characteristics:
#     results/trial_design/IND_trial_log.txt
#
# Note:
#   FORK-based parallel processing is intended for macOS and Linux.
# -------------------------------------------------------------------------



# ---- Set root directory and load files ----

project_root <- "YOURPATH/BHARP_release" # Please modify this path.
# If user entered a wrong path, stop early
if (!dir.exists(project_root)) {
  stop("ERROR: project_root does not exist. Please modify the project_root path at the top of the script.")
}

source(   file.path(project_root, "functions", "IND_model.R"))
source(   file.path(project_root, "functions", "helper_trial.R"))
source(   file.path(project_root, "scripts", "trial_design_config.R"))

IND_results_dir  <- file.path(project_root,"results", "trial_design","trial_IND")
if (!dir.exists(IND_results_dir)) {dir.create(IND_results_dir, showWarnings = FALSE, recursive = TRUE)}
options(mc.cores = 1)




log_file <- log_file <- file.path(dirname(IND_results_dir), "IND_trial_log.txt")
log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
}
  
  
log_message("====== IND Simulation setting ======")


log_message(paste0("N_sim = ", N_sim, ", N_subgrp = ", N_subgrp,
                   ", N_arm = ", N_arm, ", L = ", L))
log_message(paste0("x_fut = ", x_fut, ", x_eff = ", x_eff))

log_message(paste0("P_fut = ", paste(P_fut, collapse = ", "),
                   ", P_eff = ", paste(P_eff, collapse = ", ")))

log_message(paste0("TotalSampleSize = { " , paste(TotalSampleSize, collapse = ", ") , " }"))


log_message("Theta matrix (true cell means):")
log_message(paste(capture.output(print(Theta)), collapse = "\n"))

log_message("IND Hyperparameters:")
log_message(paste(capture.output(str(IND_hyperparam)), collapse = "\n"))


# ----Run simulated trials----
start<-Sys.time() 

cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 2L),type="FORK")

IND_results <- parallel::parLapply(cl, 1:N_sim, function(v) {
  trial_seed <- 1110000 + v
  set.seed(trial_seed)
  
  tryCatch(
    {
      out <- IND_simulate_one_trial(
        trial_index = v,
        N_subgrp = N_subgrp,
        N_arm = N_arm,
        L = L,
        TotalSampleSize = TotalSampleSize,
        Theta = Theta,
        hyperparam_list = IND_hyperparam,
        bound_fut = x_fut,
        bound_eff = x_eff,
        P_fut = P_fut,
        P_eff = P_eff,
        base_dir = IND_results_dir,
        trial_seed = trial_seed
      )
      out$trial_index <- v
      out$seed <- trial_seed
      out
    }, 
    error=function(e) {
      return(list(error=TRUE, trial_index = v, message=e$message))
    }
   )
  }
)

stopCluster(cl);

scenario_time <- as.numeric(difftime(Sys.time(), start, units = "mins"))
scenario_time <- round(scenario_time, 2)
log_message(paste0("Total runtime: ", scenario_time, " minutes"))







successful <- sapply(
  IND_results,
  function(x) is.null(x$error)
)

IND_good_results <- IND_results[successful]
IND_bad_results <- IND_results[!successful]
n_good <- length(IND_good_results)

if (n_good == 0L) {
  stop("All IND trial simulations failed.")
}

saveRDS(IND_good_results, file = file.path(IND_results_dir,"IND_good_results_list.rds"))
log_message(paste0("Failed trials: ", length(IND_bad_results), " / ", N_sim))

IND_thetamedian <- do.call(rbind, lapply(IND_good_results, "[[", "thetamedian"))
IND_thetaq750 <- do.call(rbind,lapply(IND_good_results, "[[", "thetaq750"))
IND_thetaq250 <- do.call(rbind,lapply(IND_good_results, "[[", "thetaq250"))

IND_thetapostvar <- do.call(rbind, lapply(IND_good_results, function(x) {
  apply(x$theta, 2, var)
}))



IND_efficacy <- do.call(rbind,lapply(IND_good_results, "[[", "final_efficacy"))
IND_futility <- do.call(rbind,lapply(IND_good_results, "[[", "final_futility"))
IND_cell_sample_size <- do.call(rbind,lapply(IND_good_results, "[[", "cell_sample_size"))



# ----Calculate OCs----
# Mean Absolute Error
MAE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    mean( abs( (IND_thetamedian[,x] - ThetaVec[x])  ) )  
),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Mean Abs Error:")
log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))




#Variance of Estimates
VAR<- matrix(apply(IND_thetamedian,2,var),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("VAR of posterior median:")
log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))


#Root mean squared error
RMSE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    sqrt(mean(  (IND_thetamedian[,x] - ThetaVec[x])^2 )  )
),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("RMSE:")
log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))

#Average Inter-Quantile Ranges
IQR<-matrix(colMeans(IND_thetaq750-IND_thetaq250),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Inter-Quantile Ranges:")
log_message(paste(capture.output(print(round(IQR,3))), collapse = "\n"))



PostVar<-matrix(colMeans(IND_thetapostvar),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Posterior Variance:")
log_message(paste(capture.output(print(round(PostVar,3))), collapse = "\n"))







# Expected Cell Sample Size

CellSampleSize<-matrix(colMeans(IND_cell_sample_size),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Expected Cell Sample Size:")
log_message(paste(capture.output(print(round(CellSampleSize))), collapse = "\n"))



# Global False Positive Rate


true_effective <- ThetaVec>x_eff
true_effective_mat<-matrix(rep(true_effective,each=n_good),nrow = n_good,ncol = N_arm*N_subgrp)
fp_mat=(!true_effective_mat)& IND_efficacy
GFPR<-mean(apply(fp_mat, 1, any))
log_message("Global False Positive Rate:")
log_message(GFPR)



# Generalized Power by subgrp
GP_grp<-mean(sapply(1:n_good, function(v) 
  trial_success_subgroup(IND_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by subgrp:")
log_message(GP_grp)

# Generalized Power by arm
GP_arm<-mean(sapply(1:n_good, function(v) 
  trial_success_arm(IND_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by arm:")
log_message(GP_arm)

