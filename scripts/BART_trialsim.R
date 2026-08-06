# -------------------------------------------------------------------------
# BART Adaptive Trial Design Evaluation
#
# This script evaluates the operating characteristics of a Bayesian additive
# regression tree (BART) adaptive trial design through repeated simulations.
#
# Each trial performs stage-wise recruitment, BART estimation, and interim
# or final futility and efficacy decisions.
#
# Before running:
#   - Set project_root to the local path of the cloned repository.
#   - Review the settings in scripts/trial_design_config.R.
#   - Install the BART package.
#
# Dependencies:
#   functions/BART_model.R
#   functions/helper_trial.R
#   scripts/trial_design_config.R
#
# Output:
#   Trial-level results:
#     results/trial_design/trial_BART/
#
#   Summary operating characteristics:
#     results/trial_design/BART_trial_log.txt
#
# Note:
#   FORK-based parallel processing is intended for macOS and Linux.
# -------------------------------------------------------------------------


# ---- Set root directory and load files ----

project_root <- "YOURPATH/BHARP_release" # Please modify this path.

if (!dir.exists(project_root)) {
  stop("ERROR: project_root does not exist. Please modify the project_root path at the top of the script.")
}


source(   file.path(project_root, "functions", "BART_model.R"))
source(   file.path(project_root, "functions", "helper_trial.R"))
source(   file.path(project_root, "scripts", "trial_design_config.R"))

BART_results_dir  <- file.path(project_root, "results","trial_design","trial_BART")
if (!dir.exists(BART_results_dir)) {dir.create(BART_results_dir, showWarnings = FALSE, recursive = TRUE)}





log_file <- file.path(dirname(BART_results_dir), "BART_trial_log.txt")

log_message <- function(msg) {
  cat(paste0( msg, "\n"), file = log_file, append = TRUE)
}


log_message("====== BART Simulation setting ======")


log_message(paste0("N_sim = ", N_sim, ", N_subgrp = ", N_subgrp,
                   ", N_arm = ", N_arm, ", L = ", L))
log_message(paste0("x_fut = ", x_fut, ", x_eff = ", x_eff))

log_message(paste0("P_fut = ", paste(P_fut, collapse = ", "),
                   ", P_eff = ", paste(P_eff, collapse = ", ")))

log_message(paste0("TotalSampleSize = { " , paste(TotalSampleSize, collapse = ", ") , " }"))


log_message("Theta matrix (true cell means):")
log_message(paste(capture.output(print(Theta)), collapse = "\n"))
  




# ---- Run simulated trials ----



start<-Sys.time() 
cl <- makeCluster(max(1L, parallel::detectCores() - 2L),type="FORK")

BART_results <- parLapply(cl, 1:N_sim, function(v) {
  trial_seed <- 4440000 + v
  set.seed(trial_seed)
  tryCatch({
    out <- BART_simulate_one_trial(
      trial_index = v,
      N_subgrp = N_subgrp,
      N_arm = N_arm,
      L = L,
      TotalSampleSize = TotalSampleSize,
      Theta = Theta,
      bound_fut = x_fut,
      bound_eff = x_eff,
      P_fut = P_fut,
      P_eff = P_eff,
      base_dir = BART_results_dir,
      trial_seed = trial_seed
    )
    out$trial_index <- v
    out$seed <- trial_seed
    out
  }, error=function(e) {
    return(list(error=TRUE,trial_index = v, message=e$message))
  })
})

parallel::stopCluster(cl)

scenario_time <- as.numeric(difftime(Sys.time(), start, units = "mins"))
scenario_time <- round(scenario_time, 2)
log_message(paste0("Total runtime: ", scenario_time, " minutes"))




successful <- vapply(
  BART_results,
  function(x) is.null(x$error),
  logical(1)
)

BART_good_results <- BART_results[successful]
BART_bad_results <- BART_results[!successful]
n_good <- length(BART_good_results)

saveRDS(BART_good_results, file = file.path(BART_results_dir,"BART_good_results_list.rds"))
log_message(paste0("Failed trials: ", length(BART_bad_results), " / ", N_sim))


BART_thetamedian <- do.call(rbind, lapply(BART_good_results, "[[", "thetamedian"))
BART_thetaq750 <- do.call(rbind,lapply(BART_good_results, "[[", "thetaq750"))
BART_thetaq250 <- do.call(rbind,lapply(BART_good_results, "[[", "thetaq250"))

BART_thetapostvar <- do.call(rbind, lapply(BART_good_results, function(x) {
  apply(x$theta, 2, var)
}))



BART_efficacy <- do.call(rbind,lapply(BART_good_results, "[[", "final_efficacy"))
BART_futility <- do.call(rbind,lapply(BART_good_results, "[[", "final_futility"))
BART_cell_sample_size <- do.call(rbind,lapply(BART_good_results, "[[", "cell_sample_size"))





# ----Calculate OCs----
# Mean Absolute Error
MAE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    mean( abs( (BART_thetamedian[,x] - ThetaVec[x])  ) )  
),nrow = N_arm,ncol=N_subgrp,byrow = T)


log_message("Mean Abs Error:")
log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))




# Variance of Estimates
VAR<- matrix(apply(BART_thetamedian,2,var),nrow = N_arm,ncol=N_subgrp,byrow = T)

log_message("VAR of posterior median:")
log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))


# Root mean squared error
RMSE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    sqrt(mean(  (BART_thetamedian[,x] - ThetaVec[x])^2 )  )
),nrow = N_arm,ncol=N_subgrp,byrow = T)

log_message("RMSE:")
log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))







# Average Inter-Quantile Ranges
IQR<-matrix(colMeans(BART_thetaq750-BART_thetaq250,na.rm = TRUE),nrow = N_arm,ncol=N_subgrp,byrow = T)

log_message("Inter-Quantile Ranges:")
log_message(paste(capture.output(print(round(IQR,3))), collapse = "\n"))




PostVar<-matrix(colMeans(BART_thetapostvar, na.rm = TRUE),nrow = N_arm,ncol=N_subgrp,byrow = T)


log_message("Posterior Variance:")
log_message(paste(capture.output(print(round(PostVar,3))), collapse = "\n"))








# Expected Cell Sample Size

CellSampleSize<-matrix(colMeans(BART_cell_sample_size, na.rm = TRUE),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Expected Cell Sample Size:")
log_message(paste(capture.output(print(round(CellSampleSize))), collapse = "\n"))



# Global False Positive Rate


true_effective <- ThetaVec>x_eff

true_effective_mat<-matrix(rep(true_effective,each=n_good),nrow = n_good,ncol = N_arm*N_subgrp)

fp_mat=(!true_effective_mat)& BART_efficacy

GFPR<-mean(apply(fp_mat, 1, any))
log_message("Global False Positive Rate:")
log_message(GFPR)



# Generalized Power by subgrp
GP_grp<-mean(sapply(1:n_good, function(v) 
  trial_success_subgroup(BART_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by subgrp:")
log_message(GP_grp)


# Generalized Power by arm
GP_arm<-mean(sapply(1:n_good, function(v) 
  trial_success_arm(BART_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))


log_message("Generalized Power by arm:")
log_message(GP_arm)





