# -------------------------------------------------------------------------
# BART Model and Trial Simulation Functions
#
# This file implements the Bayesian additive regression tree (BART) comparator
# used in the simulation study.
#
# Separate BART models are fitted for each treatment arm. Biomarker subgroups
# are represented using indicator variables, and posterior predictions are
# obtained for every arm-subgroup combination.
#
# Main functions:
#   - BART_analysis(): fits BART to a single dataset.
#   - BART_simulate_one_trial(): simulates and analyzes one adaptive trial.
#
# Dependencies:
#   - BART
#   - recruit_patients(), deactivate_decision(), and final_decision()
#     from functions/helper_trial.R
# -------------------------------------------------------------------------



BART_analysis<-function(Data,n_subgrp,n_arm,ndpost=8000){
    Data$k <- factor(Data$k, levels = seq_len(n_subgrp))
    k_dummies <- model.matrix(~ k - 1, data = Data)  
    colnames(k_dummies) <- paste0("k", seq_len(n_subgrp))
    Data <- cbind(Data, k_dummies)
    
    X_test <- diag(n_subgrp)
    colnames(X_test) <- paste0("k", seq_len(n_subgrp))
    theta_by_arm <- vector("list", n_arm)
    
    for (i in seq_len(n_arm)) {
    
      Dat_i <- Data[Data$i == i, , drop = FALSE]
      X_train <- as.matrix(
        Dat_i[, paste0("k", seq_len(n_subgrp)), drop = FALSE]
      )
      Y_train <- as.numeric(Dat_i$Y)
    
      bart_fit <- BART::wbart(
        x.train = X_train,
        y.train = Y_train,
        x.test  = X_test,
        ndpost  = ndpost
      )
    

      theta_by_arm[[i]] <- bart_fit$yhat.test
    }
    # Arm-major ordering:
    # Arm1Grp1, ..., Arm1GrpK, Arm2Grp1, ..., ArmIGrpK.
    theta_post <- do.call(cbind, theta_by_arm)
    
    colnames(theta_post) <- unlist(
      lapply(seq_len(n_arm), function(i) {
        paste0("Arm", i, "Grp", seq_len(n_subgrp))
      })
    )
    
    return(theta_post)
}

BART_simulate_one_trial <- function(
    trial_index, N_subgrp, N_arm,
    L,                 
    TotalSampleSize,     
    Theta,               
    bound_fut, bound_eff,          # futility / efficacy boundaries
    P_fut, P_eff,         # futility / efficacy thresholds
    base_dir,             
    trial_seed = NULL
){
  

  sim_dir <- file.path(base_dir, paste0("Sim", trial_index))
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  }
  data_dir <- file.path(sim_dir, "Data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  Status_Log <- matrix(NA, nrow=L, ncol=N_subgrp*N_arm)
  cell_actflag <- matrix(TRUE,  nrow=N_arm, ncol=N_subgrp)
  cell_final_efficacy <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  cell_final_futility <- matrix(FALSE, nrow=N_arm, ncol=N_subgrp)
  Dat <- data.frame(i=integer(0), k=integer(0), Y=numeric(0))
  

  thetamedian <- rep(NA, N_subgrp*N_arm)
  thetaq750   <- rep(NA, N_subgrp*N_arm)
  thetaq250   <- rep(NA, N_subgrp*N_arm)
  final_efficacy_vec <- rep(NA, N_subgrp*N_arm)
  final_futility_vec <- rep(NA, N_subgrp*N_arm)

  for (l in seq_len(L)) {

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
    
    if (!is.null(trial_seed)) {
      set.seed(trial_seed + l)
    }
    
    theta_l <- BART_analysis(
      Data = Dat,
      n_subgrp = N_subgrp,
      n_arm = N_arm
    )
    
    theta_l <- as.matrix(theta_l)
    if (is.null(dim(theta_l))) {
      stop("theta_l has no dim at trial ", trial_index, ", stage ", l)
    }
    
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
      

      if (all(!cell_actflag)) {

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
        
        
        
        break  
      }
      
    } 
    else {
      # l == L (final)
      write.csv(theta_l, file.path(sim_dir, paste0("theta_stage", l, ".csv")),row.names = FALSE)
      thetamedian <- apply(theta_l, 2, median)
      thetaq750   <- apply(theta_l, 2, quantile, prob=0.75)
      thetaq250   <- apply(theta_l, 2, quantile, prob=0.25)
      
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
      
      
      
    } 
    
    
    
  } 
  
  return(list(
    sample_size=TotalSampleSize[l],
    cell_sample_size=as.vector(t(cell_sample_size)),
    theta=theta_l,
    thetamedian=thetamedian,
    thetaq750=thetaq750,
    thetaq250=thetaq250,
    final_efficacy=final_efficacy_vec,
    final_futility=final_futility_vec
  ))
}

