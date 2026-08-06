# -------------------------------------------------------------------------
# BLAST Adaptive Trial Design Evaluation
#
# This script evaluates the operating characteristics of the BLAST
# adaptive trial design through repeated trial simulations.
#
# For each simulated trial, the script:
#   1. recruits patients according to the prespecified sample sizes;
#   2. performs interim and final BLAST analyses;
#   3. applies the shared futility and efficacy decision rules; and
#   4. records estimation and design operating characteristics.
#
# Before running:
#   - Set project_root to the local path of the cloned repository.
#   - Install the R packages listed under "Dependencies" below.
#   - Review the settings in scripts/trial_design_config.R.
#
# Dependencies:
#   R packages:
#     rstan
#
#   The parallel package is included with R.
#
#   Repository files:
#     functions/BLAST_model.R
#     functions/helper_cocluster.R
#     functions/helper_trial.R
#     scripts/trial_design_config.R
#
# Output:
#   Trial-level results are written to:
#
#     results/trial_design/trial_BLAST/
#
#   Summary operating characteristics are written to:
#
#     results/trial_design/BLAST_trial_log.txt
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


source(file.path(project_root, "functions", "helper_cocluster.R"))
source(file.path(project_root, "functions", "BLAST_model.R"))
source(file.path(project_root, "functions", "helper_trial.R"))
source(file.path(project_root, "scripts", "trial_design_config.R"))


BLAST_results_dir <- file.path(
  project_root,
  "results",
  "trial_design",
  "trial_BLAST"
)
if (!dir.exists(BLAST_results_dir)) {dir.create(BLAST_results_dir, showWarnings = FALSE, recursive = TRUE)}
options(mc.cores = 1)




# Creating simulation log
log_file <-  file.path(dirname(BLAST_results_dir), "BLAST_trial_log.txt")
  
log_message <- function(msg) {
    cat(paste0( msg, "\n"), file = log_file, append = TRUE)
  }
  
log_message("====== BLAST Simulation setting ======")
log_message(paste0("N_sim = ", N_sim, ", N_subgrp = ", N_subgrp,
                     ", N_arm = ", N_arm, ", L = ", L))
log_message(paste0("x_fut = ", x_fut, ", x_eff = ", x_eff))
  
log_message(paste0("P_fut = ", paste(P_fut, collapse = ", "),
                     ", P_eff = ", paste(P_eff, collapse = ", ")))
  
log_message(paste0("TotalSampleSize = { " , paste(TotalSampleSize, collapse = ", ") , " }"))

log_message("Theta matrix (true cell means):")
log_message(paste(capture.output(print(Theta)), collapse = "\n"))
  
log_message("BLAST Hyperparameters:")
log_message(paste(capture.output(str(BLAST_hyperparam)), collapse = "\n"))
  
  
  

# ---- Run simulated trials ----

start <- Sys.time() 

cl <- parallel::makeCluster(max(1L,parallel::detectCores() - 2L), type="FORK")

BLAST_results <- parallel::parLapply(
  cl, 
  1:N_sim, 
  function(v) {
    trial_seed <- 3330000 + v
    set.seed(trial_seed)
    tryCatch(
      {

        out <- BLAST_simulate_one_trial(
          trial_index = v, 
          N_subgrp = N_subgrp, 
          N_arm = N_arm,
          L = L,
          TotalSampleSize = TotalSampleSize,
          Theta = Theta,
          hyperparam_list = BLAST_hyperparam,
          bound_fut = x_fut, 
          bound_eff = x_eff, 
          P_fut = P_fut, 
          P_eff = P_eff, 
          base_dir = BLAST_results_dir,
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

parallel::stopCluster(cl)

scenario_time <- as.numeric(difftime(Sys.time(), start, units = "mins"))
scenario_time <- round(scenario_time, 2)
log_message(paste0("Total runtime: ", scenario_time, " minutes"))


successful <- sapply(
  BLAST_results,
  function(x) {
    is.null(x$error)
  }
)

BLAST_good_results <- BLAST_results[successful]
BLAST_bad_results <- BLAST_results[!successful]

n_good <- length(BLAST_good_results)
if (n_good == 0L) {
  stop("All BLAST trial simulations failed.")
}
saveRDS(BLAST_good_results, file = file.path(BLAST_results_dir,"BLAST_good_results_list.rds"))

log_message(paste0("Failed trials: ", length(BLAST_bad_results), " / ", N_sim))


 
BLAST_thetamedian <- do.call(rbind, lapply(BLAST_good_results, "[[", "thetamedian"))
BLAST_thetaq750 <- do.call(rbind,lapply(BLAST_good_results, "[[", "thetaq750"))
BLAST_thetaq250 <- do.call(rbind,lapply(BLAST_good_results, "[[", "thetaq250"))
BLAST_thetapostvar <- do.call(rbind, lapply(BLAST_good_results, function(x) { apply(x$theta, 2, var) }))



BLAST_efficacy <- do.call(rbind,lapply(BLAST_good_results, "[[", "final_efficacy"))
BLAST_futility <- do.call(rbind,lapply(BLAST_good_results, "[[", "final_futility"))
BLAST_cell_sample_size <- do.call(rbind,lapply(BLAST_good_results, "[[", "cell_sample_size"))



BLAST_q_list <- lapply(BLAST_good_results, function(x) x$q)




# Mean Absolute Error
MAE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x)
    mean( abs( (BLAST_thetamedian[,x] - ThetaVec[x])  ) )
),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Mean Abs Error:")
  log_message(paste(capture.output(print(round(MAE,3))), collapse = "\n"))




# Variance of Estimates
VAR<- matrix(apply(BLAST_thetamedian,2,var),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("VAR of posterior median:")
log_message(paste(capture.output(print(round(VAR,3))), collapse = "\n"))


# Root mean squared error
RMSE<-matrix(sapply(
  seq_len(N_arm*N_subgrp), function(x) 
    sqrt(mean(  (BLAST_thetamedian[,x] - ThetaVec[x])^2 )  )
),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("RMSE:")
log_message(paste(capture.output(print(round(RMSE,3))), collapse = "\n"))



# Average Inter-Quantile Ranges
IQR<-matrix(colMeans(BLAST_thetaq750-BLAST_thetaq250,na.rm = TRUE),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Inter-Quantile Ranges:")
log_message(paste(capture.output(print(round(IQR,3))), collapse = "\n"))

# Posterior Variance

PostVar<-matrix(colMeans(BLAST_thetapostvar,na.rm = TRUE),nrow = N_arm,ncol=N_subgrp,byrow = T)
log_message("Posterior Variance:")
log_message(paste(capture.output(print(round(PostVar,3))), collapse = "\n"))



# Expected Cell Sample Size

CellSampleSize<-matrix(colMeans(BLAST_cell_sample_size,na.rm = TRUE),nrow = N_arm,ncol=N_subgrp,byrow = T)

{
  log_message("Expected Cell Sample Size:")
  log_message(paste(capture.output(print(round(CellSampleSize))), collapse = "\n"))
}


# Global False Positive Rate

true_effective <- ThetaVec>x_eff

true_effective_mat<-matrix(rep(true_effective,each=n_good),nrow = n_good,ncol = N_arm*N_subgrp)

fp_mat=(!true_effective_mat)& BLAST_efficacy

GFPR<-mean(apply(fp_mat, 1, any))
log_message("Global False Positive Rate:")
log_message(GFPR)

# Generalized Power by subgrp
GP_grp<-mean(sapply(1:n_good, function(v) 
  trial_success_subgroup(BLAST_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by subgrp:")
log_message(GP_grp)


# Generalized Power by arm
GP_arm<-mean(sapply(1:n_good, function(v) 
  trial_success_arm(BLAST_efficacy[v,], true_effective, N_subgrp=N_subgrp, N_arm=N_arm)))
log_message("Generalized Power by arm:")
log_message(GP_arm)


# q

BLAST_q_mat <- do.call(rbind, BLAST_q_list)
colnames(BLAST_q_mat) <- paste0("Arm", seq_len(N_arm))

BLAST_q_dist_by_arm <- lapply(seq_len(N_arm), function(i) {
  round(prop.table(table(BLAST_q_mat[, i])), 3)
})

names(BLAST_q_dist_by_arm) <- paste0("Arm", seq_len(N_arm))

log_message("Distribution of selected q for BLAST by arm:")
log_message(paste(capture.output(print(BLAST_q_dist_by_arm)), collapse = "\n"))


