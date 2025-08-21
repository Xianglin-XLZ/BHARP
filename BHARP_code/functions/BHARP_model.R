# File: BHARP_model.R
# Purpose:
#   This script implements the core analysis and trial simulation components of
#   the BHARP (Bayesian Hierarchical Adjustable Random Partition) method.
#
# Contents:
#   - ReadCombineMCMC: Read and combine MCMC samples from multiple chains/files.
#   - CoClusterProb / ArmListCoCluster: Compute co-cluster probabilities from latent partitions.
#   - Bharp_analysis: Perform BHARP MCMC analysis at a given interim/final stage.
#   - simulate_one_trial: Simulate a complete adaptive trial using the BHARP model,
#                         with interim analyses, early stopping for efficacy/futility,
#                         and posterior summarization.
#
# Dependencies:
#   `Bharp()` is defined externally in Rcpp.
#   Relies on external functions:
#     - recruit_patients(): for data generation at each analysis stage
#     - deactivate_decision(), final_decision(): for adaptive decision-making




# to read and combine MCMC sample 
ReadCombineMCMC<- function(prefix,header=FALSE) {
  #prefix can include path/param_dataid-
  dir_path <- dirname(prefix)  
  file_prefix <- basename(prefix)  
  if (dir_path == ".") {
    dir_path <- getwd()
  }
  file_pattern <- paste0("^", file_prefix, ".*\\.csv$")
  csv_files <- list.files(path = dir_path, pattern = file_pattern, full.names = TRUE)
  
  if (length(csv_files) > 0) {
    combined_data <- do.call(rbind, lapply(csv_files, function(file) {
      read.csv(file, header = header)
    }))
    return(combined_data)
  } else {
    warning(paste("No files found in", dir_path, "matching pattern:", file_pattern))
    return(NULL)
  }
}


# to calculate co-cluster probability for one arm
CoClusterProb<-function(iz){
  s <- ncol(iz)     # Narm
  Prob_matrix <- matrix(NA, nrow = s, ncol = s)   # Matrix of co-cluster probability
  
  for (i in 1:(s-1)) {
    for (j in (i+1):s) {
      Prob_matrix[i, j] <- mean(iz[, i] == iz[, j]) 
    }
  }
  diag(Prob_matrix)<-1
  return(round(Prob_matrix,4))
}

# to read MCMC output and calculate pairwise coCluster probability for all arms
ArmListCoCluster<-function(combined_iz,Nsubgrp,Narm){
  
  list_iz<- list()  
  for (i in 1:Narm){
    col_start <- (i - 1) * Nsubgrp + 1  
    col_end <- col_start + Nsubgrp - 1  
    list_iz[[i]] <- combined_iz[, col_start:col_end, drop = FALSE]
  }    #split iz by arm
  
  probability_list <- lapply(list_iz, CoClusterProb)
  return(probability_list)
}





# conduct analysis MCMC (sample =2000 thinning=3 burn=1000) *4chains
# return posterior sample of theta and iz for decisions
Bharp_analysis<-function(l,dataset,n_subgrp,n_arm,
                         hyperparam_list,results_dir){
  # l: analysis stage
  # dataset: accumulated dataset when conducting the analysis 
  # hyperparam_list: contains hyperparameters of the prior
  # results_dir: the folder saving MCMC output
  # bound, prob: numbers saving the boundaries
  
  Bharp(dataset$i, dataset$k, dataset$Y, ngrp=n_subgrp, narm=n_arm,
        c_beta=hyperparam_list$c_beta,       p_beta=hyperparam_list$p_beta,
        a_cell=hyperparam_list$a_cell,       b_cell=hyperparam_list$b_cell,
        a_between=hyperparam_list$a_between, b_between=hyperparam_list$b_between,
        a_within=hyperparam_list$a_within,   b_within=hyperparam_list$b_within,
        results_dir=results_dir,
        nsamp=2000,thn=3,burn=1000,dataid=l,nchain=4)
  
  #get parameter posterior by read combineMCMC
  betaPost<-ReadCombineMCMC(paste0(results_dir,"/beta_",l,"-"))
  DeltaPost<-ReadCombineMCMC(paste0(results_dir,"/Delta_",l,"-"))
  izPost<-as.matrix(ReadCombineMCMC(paste0(results_dir,"/iz_",l,"-")))
  
  thetaPost<-as.matrix(DeltaPost+betaPost[,rep(1:n_arm,each=n_subgrp)])
  
  return(list(theta_posterior=thetaPost,iz_posterior=izPost,beta_posterior=betaPost,Delta_posterior=DeltaPost))
  
}  




# function to simulate one trial
simulate_one_trial <- function(
    trial_index, N_subgrp, N_arm,
    L,                   # number of analyses
    TotalSampleSize,     # A vector containing total sample size at each analysis
    Theta,               # Real Theta
    hyperparam_list,  
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
    
    # get theta and iz
    analysis_result_l <- Bharp_analysis(l,
                                        dataset=Dat,
                                        n_subgrp=N_subgrp,
                                        n_arm=N_arm,
                                        hyperparam=hyperparam_list, 
                                        results_dir    )
    theta_l <- analysis_result_l$theta_posterior
    iz_l    <- analysis_result_l$iz_posterior
    
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
        
        write.csv(iz_l, file.path(sim_dir, paste0("iz_stage", l, ".csv")),row.names = FALSE)
        iz_cocluster_prob<-ArmListCoCluster(iz_l, Nsubgrp = N_subgrp, Narm=N_arm)
        
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
      #save iz
      write.csv(iz_l, file.path(sim_dir, paste0("iz_stage", l, ".csv")),row.names = FALSE)
      iz_cocluster_prob<-ArmListCoCluster(iz_l, Nsubgrp = N_subgrp, Narm=N_arm)
      
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
    iz_cocluster_prob=iz_cocluster_prob,
    final_efficacy=final_efficacy_vec,
    final_futility=final_futility_vec,
    avg_effect_nonfutile=avg_effect_nonfutile,
    avg_effect_effective=avg_effect_effective
  ))
}
