# -------------------------------------------------------------------------
# BLAST Model and Trial Simulation Functions
#
# This file implements BLAST as a finite mixture model for subgroup-level
# treatment effects. The number of mixture components is selected from
# q = 1, 2, 3 using the Deviance Information Criterion (DIC).
#
# Functions and objects:
#   BLAST_code:
#     Stan specification of the BLAST finite mixture model.
#
#   BLAST_stan_model:
#     Compiled Stan model.
#
#   dic():
#     Calculates DIC from a fitted BLAST model.
#
#   BLAST_analysis():
#     Fits the candidate single-arm models and selects q using DIC.
#
#   BLAST_analysis_multiarm():
#     Applies BLAST independently to each treatment arm.
#
#   BLAST_simulate_one_trial():
#     Simulates one adaptive trial using the BLAST analyses.
#
# Dependencies:
#   R package:
#     rstan
#
#   Repository functions:
#     ArmListCoCluster() is defined in functions/helper_cocluster.R.
#
#     recruit_patients(), deactivate_decision(), and final_decision() are
#     defined in functions/helper_trial.R.
# -------------------------------------------------------------------------



# Stan model


BLAST_code <- "

data {
  int<lower=0> N;                  
  int<lower=1> Nsubgrp;           
  vector[N] Y;                     
  int<lower=1,upper=Nsubgrp> kIndex[N];  
  int<lower=1,upper=Nsubgrp> q;   
  
  real a_cell;
  real b_cell;
  real a_between;
  real b_between;
  real a_within;
  real b_within;
}

parameters {
  real<lower=0> varsigma;     
  real<lower=0> tau;          
  
  simplex[q] w;               
  vector<lower=0>[q] sigma;   
  vector[q] mu;              
  
  vector[Nsubgrp] theta;      
}

model {
  // Priors (aligned with BHARP)
  varsigma ~ gamma(a_cell, b_cell);          
  tau      ~ gamma(a_between, b_between);     
  mu       ~ normal(0, inv_sqrt(tau));        
  sigma    ~ inv_gamma(a_within, b_within);   
  w        ~ dirichlet(rep_vector(1.0, q));   
  
   // Finite mixture prior for subgroup-level treatment effects
  for (k in 1:Nsubgrp) {
    vector[q] lp_mix;
    for (t in 1:q)
      lp_mix[t] = log(w[t]) + normal_lpdf(theta[k] | mu[t], sqrt(sigma[t]));
    target += log_sum_exp(lp_mix);
  }
  
  // Outcome likelihood
  Y ~ normal(theta[kIndex], inv_sqrt(varsigma));
}

generated quantities {
  vector[N] log_lik;
  int z[Nsubgrp];                  
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(Y[n] | theta[kIndex[n]], inv_sqrt(varsigma));
  }

  // Draw the posterior component allocation for each subgroup
  for (k in 1:Nsubgrp) {
    vector[q] lp_mix;
    for (t in 1:q) {
      lp_mix[t] = log(w[t]) + normal_lpdf(theta[k] | mu[t], sqrt(sigma[t]));
    }
    lp_mix = lp_mix - log_sum_exp(lp_mix);
    z[k] = categorical_rng(softmax(lp_mix));
  }
}


"

BLAST_stan_model <- rstan::stan_model(model_code = BLAST_code)

# DIC calculation
dic <- function(fit, standata) {
  samples <- rstan::extract(
    fit,
    pars = c("log_lik", "theta", "varsigma")
  )
  log_lik <- samples$log_lik
  theta_draws <- samples$theta
  varsigma_draws <- samples$varsigma
  
  kIndex <- as.integer(standata$kIndex)
  Y <- as.numeric(standata$Y)
  # basic checks
  if (length(dim(log_lik)) != 2) {
    stop("log_lik should be a draws x observations matrix.")
  }
  if (ncol(log_lik) != length(Y)) {
    stop("ncol(log_lik) does not match length(Y).")
  }
  if (length(dim(theta_draws)) != 2) {
    stop("theta should be a draws x subgroups matrix. If this is multi-arm joint model, arm-specific indexing is needed.")
  }
  if (length(kIndex) != length(Y)) {
    stop("length(kIndex) does not match length(Y).")
  }
  if (max(kIndex) > ncol(theta_draws) || min(kIndex) < 1) {
    stop("kIndex is out of range for theta.")
  }
  dev_bar <- -2 * mean(rowSums(log_lik))
  theta_mean <- colMeans(theta_draws)
  varsigma_mean <- mean(varsigma_draws)
  mu_y <- theta_mean[kIndex]
  sd_y <- 1 / sqrt(varsigma_mean)
  loglik_mean <- dnorm(
    Y,
    mean = mu_y,
    sd = sd_y,
    log = TRUE
  )
  if (
    !is.finite(dev_bar) ||
    !is.finite(varsigma_mean) ||
    !is.finite(sd_y) ||
    any(!is.finite(theta_mean)) ||
    any(!is.finite(mu_y)) ||
    any(!is.finite(loglik_mean))
  ) {
    dev <- NA_real_
    DIC <- NA_real_
  } else {
    dev <- -2 * sum(loglik_mean)
    DIC <- 2 * dev_bar - dev
  }
  
  data.frame(
    D_bar = dev_bar,
    D_hat = dev,
    DIC = DIC
  )
}




BLAST_analysis <- function(
    Data,
    n_subgrp,
    hyperparam,
    seed=NULL
) {
  rstan::rstan_options(auto_write = TRUE)
  standata <- list(
    N       = length(Data$Y),
    Nsubgrp = n_subgrp,
    Y       = Data$Y,
    kIndex  = Data$k,
    a_cell    = hyperparam$a_cell,
    b_cell    = hyperparam$b_cell,
    a_between = hyperparam$a_between,
    b_between = hyperparam$b_between,
    a_within  = hyperparam$a_within,
    b_within  = hyperparam$b_within
  )
  
  fit_list <- lapply(1:3, function(q_now) {
    seed_q <- if (is.null(seed)) NULL else seed + q_now
    rstan::sampling(
      BLAST_stan_model,
      data   = c(standata, list(q = q_now)),
      iter   = 3000,
      chains = 4, 
      warmup=1000, 
      seed=seed_q,
      cores =1, 
      refresh=0
    )
  }
)
  
  dic_tbl <- do.call(
    rbind,
    Map(function(f, q_now) {
      cbind(q = q_now, dic(f, standata))
    }, 
    fit_list, 
    1:3)
  )
  
  dic_tbl <- as.data.frame(dic_tbl)
  dic_tbl$q   <- as.integer(dic_tbl$q)
  dic_tbl$DIC <- as.numeric(dic_tbl$DIC)
  valid_rows <- which(is.finite(dic_tbl$DIC))
  
  if (length(valid_rows) == 0) {
    stop(
      paste0(
        "No valid DIC values for q = 1, 2, 3. DIC table:\n",
        paste(capture.output(print(dic_tbl)), collapse = "\n")
      )
    )
  }
  
  best_row <- valid_rows[which.min(dic_tbl$DIC[valid_rows])]
  best_q   <- dic_tbl$q[best_row]
  
  if (length(best_q) != 1 || is.na(best_q) || !(best_q %in% seq_along(fit_list))) {
    stop(
      paste0(
        "Invalid selected q: ", paste(best_q, collapse = ", "),
        "\nDIC table:\n",
        paste(capture.output(print(dic_tbl)), collapse = "\n")
      )
    )
  }
  
  best_fit <- fit_list[[best_q]]
  samples  <- rstan::extract(best_fit)
  
  theta_draws <- matrix(samples$theta, ncol = n_subgrp)
  z_draws     <- matrix(samples$z,     ncol = n_subgrp)
  
  if (ncol(theta_draws) != n_subgrp) {
    stop("BLAST theta has wrong shape.")
  }
  if (ncol(z_draws) != n_subgrp) {
    stop("BLAST z has wrong shape.")
  }

  iz_cocluster_prob <- ArmListCoCluster(z_draws, n_subgrp, 1)
  
  list(
    q                    = best_q,
    theta_posterior      = theta_draws,
    iz_posterior         = z_draws,
    iz_cocluster_prob    = iz_cocluster_prob,
    dic_table            = dic_tbl
  )
}

# Fit the single-arm BLAST model independently to each treatment arm.
BLAST_analysis_multiarm <- function(
    Data, 
    n_subgrp, 
    n_arm, 
    hyperparam,
    seed=NULL
) {
  res_list <- lapply(seq_len(n_arm), function(i_now) {
    arm_seed <- if (is.null(seed)) NULL else seed + 1000L * i_now
    Data_i <- Data[Data$i == i_now, c("Y", "k"), drop = FALSE]
    Data_i$k <- as.integer(Data_i$k)
    BLAST_analysis(
      Data       = Data_i,
      n_subgrp   = n_subgrp,
      hyperparam = hyperparam,
      seed       = arm_seed
    )
  }
 )
  
  theta_all <- do.call(
    cbind,
    lapply(res_list, function(x) {
      theta_i <- matrix(x$theta_posterior, ncol = n_subgrp)
      if (ncol(theta_i) != n_subgrp) {
        stop("Arm-level BLAST theta has wrong shape.")
      }
      
      theta_i
    })
  )
  
  z_all <- do.call(
    cbind,
    lapply(res_list, function(x) {
      z_i <- matrix(x$iz_posterior, ncol = n_subgrp)
      
      if (ncol(z_i) != n_subgrp) {
        stop("Arm-level BLAST z has wrong shape.")
      }
    
      z_i
    })
  )

  
  dic_table_all <- do.call(
    rbind,
    Map(function(x, i_now) {
      cbind(arm = i_now, x$dic_table)
    }, res_list, seq_len(n_arm))
  )
  
  list(
    q                 = sapply(res_list, function(x) x$q),
    theta_posterior   = theta_all,
    iz_posterior      = z_all,
    iz_cocluster_prob = lapply(res_list, function(x) x$iz_cocluster_prob),
    dic_table         = dic_table_all
  )
}



# function to simulate one trial
BLAST_simulate_one_trial <- function(
    trial_index,
    N_subgrp, 
    N_arm,
    L,                   
    TotalSampleSize,     
    Theta,               
    hyperparam_list,  
    bound_fut, bound_eff,          # futility / efficacy boundaries
    P_fut, P_eff,          # futility / efficacy thresholds
    base_dir,    
    trial_seed=NULL
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
  
  q_terminal <- rep(NA_integer_, N_arm)
  dic_table_terminal <- NULL
  

  for (l in seq_len(L)) {

    Dat <- recruit_patients(
      Dat,
      cell_active = cell_actflag,
      targetN = TotalSampleSize[l],
      n_subgrp = N_subgrp,
      n_arm = N_arm,
      Theta = Theta,
      method = "balanced"
    )
    
    write.csv(Dat, file.path(data_dir, paste0("Data_stage", l, ".csv")))
    stage_seed <- if (is.null(trial_seed)) NULL else trial_seed + 100L * l
    analysis_result_l <- BLAST_analysis_multiarm(
      Dat,
      n_subgrp = N_subgrp,
      n_arm = N_arm,
      hyperparam = hyperparam_list,
      seed = stage_seed
    )
    
    theta_l <- analysis_result_l$theta_posterior
    iz_l    <- analysis_result_l$iz_posterior
    q_l     <- analysis_result_l$q
    theta_l <- matrix(theta_l, ncol = N_arm * N_subgrp)
    
    if (is.null(dim(theta_l)) || ncol(theta_l) != N_arm * N_subgrp) {
      stop("theta_l has wrong shape at trial ", trial_index, ", stage ", l, ". dim = ", paste(dim(theta_l), collapse = " x ") )
    }
    
    q_terminal <- q_l
    dic_table_terminal <- analysis_result_l$dic_table
    
    write.csv(
      data.frame(arm = seq_len(N_arm), q = q_l),
      file.path(sim_dir, paste0("BLAST_q_stage", l, ".csv")),
      row.names = FALSE
    )
    

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
        
        write.csv(iz_l, file.path(sim_dir, paste0("iz_stage", l, ".csv")),row.names = FALSE)
    
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
      #save iz
      write.csv(iz_l, file.path(sim_dir, paste0("iz_stage", l, ".csv")),row.names = FALSE)
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
  
  # return information for one simulation
  return(list(
    sample_size=TotalSampleSize[l],
    cell_sample_size=as.vector(t(cell_sample_size)),
    theta=theta_l,
    thetamedian=thetamedian,
    thetaq750=thetaq750,
    thetaq250=thetaq250,
    q = q_terminal,
    dic_table = dic_table_terminal,
    final_efficacy=final_efficacy_vec,
    final_futility=final_futility_vec
  ))
}





