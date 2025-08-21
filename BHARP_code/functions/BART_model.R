# File: BART_model.R
# Purpose:
#   This script defines the BART-based (Bayesian Additive Regression Trees) model
#   used as a nonparametric baseline in the evaluation of treatment effect heterogeneity.
#   It includes functions for posterior estimation using the `BART` package and
#   simulation of adaptive clinical trials with subgroup-specific effects.
#
# Contents:
#   - BART_analysis(): Fit separate BART models for each treatment arm, with subgroup indicators
#                      as predictors. Returns posterior samples of subgroup-specific effects.
#   - BART_simulate_one_trial(): Simulate an adaptive multi-stage trial using BART-based
#                                effect estimation and decision rules for early stopping.
#
# Dependencies:
#   - Requires the `BART` package (e.g., `wbart()` function).
#   - Assumes input data frame contains columns: `i` (arm index), `k` (subgroup index), and `Y`.
#   - Uses external functions: 
#       - recruit_patients()            (for generating new patients by subgroup)
#       - deactivate_decision(), final_decision() (for adaptive interim and final decisions)


#----BART----


BART_analysis<-function(Data,n_subgrp,n_arm){
  
  if (n_subgrp > 1) {
    Data$k <- factor(Data$k)
    k_dummies <- model.matrix(~ k - 1, data = Data)  
    Data <- cbind(Data, k_dummies)
  } else {
    Data$k1 <- 1                                    
  }
  
  one_BART_list <- vector("list", n_subgrp * n_arm)
  
  idx <- 1
  for (i in seq_len(n_arm)) {
    
    Dat_i <- Data[Data$i == i, ]
    X_i   <- as.matrix(Dat_i[, paste0("k", 1:n_subgrp), drop = FALSE])
    Y_i   <- Dat_i$Y
    
    bart_mod   <- wbart(x.train = X_i,
                        y.train = Y_i,
                        ndpost   = 2000)
    
    yhat_mat   <- bart_mod$yhat.train   
    for (k in seq_len(n_subgrp)) {
      idx_k <- which(Dat_i$k == k)
      cell_mean_draws <- rowMeans(yhat_mat[, idx_k, drop = FALSE])
      one_BART_list[[idx]] <- cell_mean_draws
      idx <- idx + 1
    }
  }
  
  theta_post <- do.call(cbind, one_BART_list)          # ndpost × (Arm×Grp)
  theta_post
}

BART_simulate_one_trial <- function(
    trial_index, N_subgrp, N_arm,
    L,                   # number of analyses
    TotalSampleSize,     # A vector containing total sample size at each analysis
    Theta,               # Real Theta
    #hyperparam_list,  
    bound_fut, bound_eff,          # futility / efficacy boundaries
    P_fut, P_eff,          # futility / efficacy thresholds
    base_dir = getwd()     # base directory to save all simulation results
){
  
  # create folder 
  sim_dir <- file.path(base_dir, paste0("Sim", trial_index))
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  }
  results_dir <- file.path(sim_dir, "Results")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  }
  data_dir <- file.path(sim_dir, "Data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # initialize containers
  Status_Log <- matrix(NA, nrow=L, ncol=N_subgrp*N_arm)
  cell_actflag <- matrix(TRUE,  nrow=N_arm, ncol=N_subgrp)
  cell_final_efficacy <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  cell_final_futility <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  Dat <- data.frame(i=integer(0), k=integer(0), Y=numeric(0))
  
  # the vectors to return
  thetamedian <- rep(NA, N_subgrp*N_arm)
  thetaq750   <- rep(NA, N_subgrp*N_arm)
  thetaq250   <- rep(NA, N_subgrp*N_arm)
  final_efficacy_vec <- rep(NA, N_subgrp*N_arm)
  final_futility_vec <- rep(NA, N_subgrp*N_arm)
  
  # analyze by stage
  for (l in seq_len(L)) {
    # recruit and save data
    Dat <- recruit_patients(
      dataset = Dat,
      cell_active = cell_actflag,
      targetN = TotalSampleSize[l],
      n_subgrp = N_subgrp,
      n_arm = N_arm,
      Theta = Theta,
      method = "balanced"
    )
    write.csv(Dat, file.path(data_dir, paste0("Data_stage", l, ".csv")))
    
    
    theta_l <- BART_analysis(Data=Dat,
                             n_subgrp=N_subgrp,
                             n_arm=N_arm)
    
    
    # interim  vs final
    if (l < L) {
      # interim
      decision_l <- deactivate_decision(
        theta_posterior = theta_l,
        cell_active = cell_actflag,
        cell_final_efficacy = cell_final_efficacy,
        cell_final_futility = cell_final_futility,
        futility_bound = bound_fut,
        efficacy_bound = bound_eff,
        futility_threshold = P_fut[l],
        efficacy_threshold = P_eff[l]
      )
      cell_actflag        <- decision_l$Cell_Active_Updated
      cell_final_efficacy <- decision_l$Cell_Final_Efficacy_Updated
      cell_final_futility <- decision_l$Cell_Final_Futility_Updated
      Status_Log[l, ]     <- as.vector(t(decision_l$Status))
      
      # terminate trial if no active cell
      if (all(!cell_actflag)) {
        # save interim results as final
        write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
        
        
        thetamedian <- apply(theta_l, 2, median)
        thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
        thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
        
        
        final_efficacy_vec <- as.vector(t(cell_final_efficacy))
        write.csv(cell_final_efficacy, file.path(sim_dir, "final_efficacy.csv"),row.names = FALSE)
        
        final_futility_vec <- as.vector(t(cell_final_futility))
        write.csv(cell_final_futility, file.path(sim_dir, "final_futility.csv"),row.names = FALSE)
        write.csv(Status_Log, file.path(sim_dir, "status_log.csv"),row.names = FALSE)
        
        
        
        thetamedian_mat<-matrix(thetamedian,nrow = N_arm,ncol = N_subgrp,byrow=TRUE)
        
        cell_sample_size=as.matrix(table(Dat$i, Dat$k))
        #calculate nonfutile cells' weighted average estimation
        avg_effect_nonfutile <- apply(matrix(seq_len(N_arm)), 1, function(i) {
          nf_idx <- which(!cell_final_futility[i, ])
          if (length(nf_idx)==0) {
            NA
          } else {
            w <- cell_sample_size[i, nf_idx]
            sum(thetamedian_mat[i, nf_idx] * w) / sum(w)
          }
        })
        
        
        #calculate effective cells' average estimation
        avg_effect_effective <- apply(matrix(seq_len(N_arm)), 1, function(i) {
          eff_idx <- which(cell_final_efficacy[i, ])
          if (length(eff_idx)==0) {
            NA
          } else {
            w <- cell_sample_size[i, eff_idx]
            sum(thetamedian_mat[i, eff_idx] * w) / sum(w)
          }
        })
        
        
        
        break  # terminate for l in 1:L 
      }
      
    } 
    else {
      # l == L (final)
      # save theta
      write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
      thetamedian <- apply(theta_l, 2, median)
      thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
      thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
      
      # draw final conclusion for active cells
      decision_l <- final_decision(
        theta_posterior=theta_l,
        cell_active=cell_actflag,
        cell_final_efficacy=cell_final_efficacy,
        cell_final_futility=cell_final_futility,
        futility_bound=bound_fut,
        efficacy_bound=bound_eff,
        futility_threshold=P_fut[l],
        efficacy_threshold=P_eff[l]
      )
      
      cell_final_efficacy <- decision_l$Cell_Final_Efficacy_Updated
      cell_final_futility <- decision_l$Cell_Final_Futility_Updated
      Status_Log[l, ]     <- as.vector(t(decision_l$Status))
      
      final_efficacy_vec <- as.vector(t(cell_final_efficacy))
      write.csv(cell_final_efficacy, file.path(sim_dir, "final_efficacy.csv"),row.names = FALSE)
      
      final_futility_vec <- as.vector(t(cell_final_futility))
      write.csv(cell_final_futility, file.path(sim_dir, "final_futility.csv"),row.names = FALSE)
      write.csv(Status_Log, file.path(sim_dir, "status_log.csv"),row.names = FALSE)
      
      
      thetamedian_mat<-matrix(thetamedian,nrow = N_arm,ncol = N_subgrp,byrow=TRUE)
      cell_sample_size=as.matrix(table(Dat$i, Dat$k))
      #calculate nonfutile cells' weighted average estimation
      avg_effect_nonfutile <- apply(matrix(seq_len(N_arm)), 1, function(i) {
        nf_idx <- which(!cell_final_futility[i, ])
        if (length(nf_idx)==0) {
          NA
        } else {
          w <- cell_sample_size[i, nf_idx]
          sum(thetamedian_mat[i, nf_idx] * w) / sum(w)
        }
      })
      
      
      #calculate effective cells' average estimation
      avg_effect_effective <- apply(matrix(seq_len(N_arm)), 1, function(i) {
        eff_idx <- which(cell_final_efficacy[i, ])
        if (length(eff_idx)==0) {
          NA
        } else {
          w <- cell_sample_size[i, eff_idx]
          sum(thetamedian_mat[i, eff_idx] * w) / sum(w)
        }
      })
      
      
      
    } # end if (l < L) else
    
    
    
  } # end for l in seq_len(L)
  
  # return information for one simulation
  return(list(
    sample_size=TotalSampleSize[l],
    cell_sample_size=as.vector(t(cell_sample_size)),
    theta=theta_l,
    thetamedian=thetamedian,
    thetaq750=thetaq750,
    thetaq250=thetaq250,
    final_efficacy=final_efficacy_vec,
    final_futility=final_futility_vec,
    avg_effect_nonfutile=avg_effect_nonfutile,
    avg_effect_effective=avg_effect_effective
  ))
}

