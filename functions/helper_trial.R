# -------------------------------------------------------------------------
# Shared Helpers for Adaptive Trial Simulations
#
# This file defines patient recruitment, interim/final decision, and
# generalized power functions shared by all trial simulation methods.
#
# Functions:
#   random_balanced_partition():
#     Randomly divides an integer into approximately equal parts.
#
#   random_proportional_partition():
#     Divides an integer according to prespecified allocation ratios.
#
#   convert2cell():
#     Converts scalar, arm-level, or subgroup-level decision boundaries
#     to cell-level vectors.
#
#   recruit_patients():
#     Generates patient outcomes and updates the accumulated trial data.
#
#   deactivate_decision():
#     Applies interim futility and efficacy stopping rules.
#
#   final_decision():
#     Applies the final futility and efficacy decision rules.
#
#   trial_success_subgroup(), trial_success_arm():
#     Evaluate generalized power at the subgroup and arm levels.
# -------------------------------------------------------------------------




# Randomly divide an integer into approximately equal parts.
random_balanced_partition <- function(number, parts) {
  if (parts > number) {
    partition <- rep(0, parts)
    indices <- sample(1:parts, number)
    partition[indices] <- partition[indices] + 1
    
    return(partition)
  }
  
  base_size <- number %/% parts
  remainder <- number %% parts
  partition <- rep(base_size, parts)
  if (remainder > 0) {
    indices <- sample(1:parts, remainder)
    partition[indices] <- partition[indices] + 1
  }
  
  return(partition)
}

# Divide an integer according to the supplied allocation ratios.
random_proportional_partition <- function(number, ratio) {
  exact <- ratio / sum(ratio) * number
  floor_part <- floor(exact)
  remainder  <- number - sum(floor_part)
  frac_part <- exact - floor_part
  idx_order <- order(frac_part, decreasing = TRUE)
  x <- floor_part
  if (remainder > 0) {
    for (r in seq_len(remainder)) {
      x[idx_order[r]] <- x[idx_order[r]] + 1
    }
  }
  
  return(x)
}



# Convert scalar, arm-level, or subgroup-level boundaries to a cell-level
# vector in row-major order.
convert2cell<-function(number, Nrow, Ncol){
  if(length(number)==Nrow*Ncol|length(number)==1){
    return(number)
  }else if(length(number)==Nrow){
    return(rep(number,each=Ncol))
  }else if(length(number)==Ncol){
    return(rep(number,times=Nrow))
  }else {
    message("Dimension mismatch")
    return(number)
  }
}

# ---- Patient recruitment ----

# Recruit patients until the accumulated sample size reaches targetN.
#
# Arguments:
#   dataset:
#     Accumulated trial data. It may be empty at the start of a trial.
#
#   cell_active:
#     Logical matrix indicating which arm-subgroup cells are recruiting.
#
#   targetN:
#     Target cumulative sample size for the current analysis.
#
#   Theta:
#     Matrix of true cell-specific outcome means.
#
#   method:
#     Allocation method for newly recruited patients.
#
#   ratio:
#     Cell allocation ratios when method = "proportion".
recruit_patients <- function(
  dataset, 
  cell_active, 
  targetN,
  n_subgrp, 
  n_arm, 
  Theta, 
  method = c("balanced", "proportion"),
  ratio = NULL
) {
  method <- match.arg(method)
  
  # Calculate the number of patients required at the current stage.
  currentN <- nrow(dataset)
  if (is.null(currentN)) currentN <- 0  
  newN <- targetN - currentN
  if (newN <= 0) {
    message("No new participants needed. Current N >= targetN.")
    return(dataset)
  }
  
  # Identify the cells that remain open for recruitment.
  if (!all(dim(cell_active) == c(n_arm, n_subgrp))) {
    stop("cell_active should have dim: n_arm , n_subgrp.")
  }
  active_cells <- which(cell_active, arr.ind = TRUE) 
  n_active <- nrow(active_cells)  
  
  # Allocate newly recruited patients across active cells.
  if (method == "balanced") {
    new_partitions <- random_balanced_partition(newN, parts = n_active)
  } else if (method == "proportion") {
    if (is.null(ratio)) {
      ratio <- rep(1, n_active)
    } else {
      if (length(ratio) != n_active) {
        stop("length ratio !=  n_active!")
      }
    }
    new_partitions <- random_proportional_partition(newN, ratio)
  } else {
    stop("method = 'balanced' or 'proportion'.")
  }
  
  # Generate outcomes for newly recruited patients.
  new_recruits_list <- lapply(seq_len(n_active), function(idx) {
    i_val <- active_cells[idx, 1]
    k_val <- active_cells[idx, 2]
    count <- new_partitions[idx]
    
    if (count > 0) {
      data.frame(i = rep(i_val, count),
                 k = rep(k_val, count),
                 Y=rnorm(count, Theta[i_val,k_val], sd=1))
    } else {
      data.frame(i = integer(0), k = integer(0), Y=numeric(0))
    }
  })
  
  new_recruits <- do.call(rbind, new_recruits_list)
  updated_dataset <- rbind(dataset, new_recruits)
  
  return(updated_dataset)
}












# Apply interim futility and efficacy rules.
#
# Cells satisfying the futility rule are closed first. Among the remaining
# active cells, an arm is closed for efficacy only when all of its active
# cells satisfy the efficacy rule.

deactivate_decision<-function(
    theta_posterior, 
    cell_active, 
    cell_final_efficacy, 
    cell_final_futility,
    efficacy_bound, 
    futility_bound,
    efficacy_threshold, 
    futility_threshold
){
  Nrow<-nrow(cell_active)
  Ncol<-ncol(cell_active)
  Status<-matrix(NA, nrow=Nrow, ncol=Ncol)
  
  #record cells that were closed in previous analysis
  Status[!cell_active] <- "Closed" 
  cell_active_vec<-as.vector(t(cell_active)) 
  
  if(length(futility_bound)==1){
    P_futility <- apply(theta_posterior<futility_bound, 2, mean)
  } else{
    futility_bound_cell <- convert2cell(futility_bound, Nrow=Nrow, Ncol = Ncol)
    futility_bound_cell_mat <- t(replicate(nrow(theta_posterior), futility_bound_cell))
    
    if(!all(dim(futility_bound_cell_mat) == dim(theta_posterior))){ 
      message("dimension mismatch for futility bound") 
    }
    P_futility < -apply(theta_posterior<futility_bound_cell_mat, 2, mean)
  }
  
  # Close active cells that satisfy the futility rule.
  Futility_Decision_vec<-rep(NA,Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Futility_Decision_vec[x] <- P_futility[x]>futility_threshold
    if(isTRUE(Futility_Decision_vec[x])){
      cell_active_vec[x]<-FALSE 
    } 
  }
  
  Futility_Decision<-matrix(Futility_Decision_vec, nrow=Nrow, ncol=Ncol, byrow=TRUE)
  Status[Futility_Decision] <- "Close Futility"     
  cell_final_futility[Futility_Decision] <- TRUE
  
  
  
  # Calculate posterior efficacy probabilities.
  if(length(efficacy_bound)==1){
    P_efficacy<-apply(theta_posterior>efficacy_bound, 2, mean )
  } else{
    efficacy_bound_cell<-convert2cell(efficacy_bound, Nrow = Nrow, Ncol=Ncol)
    efficacy_bound_cell_mat<-t(replicate(nrow(theta_posterior), efficacy_bound_cell))
    if(!all(dim(efficacy_bound_cell_mat) == dim(theta_posterior))) {
      message("dimension mismatch for efficacy bound")
    }
    P_efficacy<-apply(theta_posterior>efficacy_bound_cell_mat, 2, mean )
  }
  
  # Evaluate efficacy among cells that remain active.
  Efficacy_Decision_vec<-rep(NA, Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Efficacy_Decision_vec[x] <- P_efficacy[x]>efficacy_threshold
   
  }
  Efficacy_Decision<-matrix(Efficacy_Decision_vec, nrow=Nrow, ncol=Ncol, byrow=TRUE)
  Status[Efficacy_Decision]<-"Active Efficacy" 
  Cell_Active_Updated<- matrix(cell_active_vec, nrow=Nrow, ncol=Ncol, byrow=TRUE)
  
  # Close an arm for efficacy when all of its remaining active cells
  # satisfy the efficacy rule.
  for (i in which(rowSums(Cell_Active_Updated) > 0) ){
    active_grp_i<-which(Cell_Active_Updated[i,]) 
    if (all( Efficacy_Decision[i, active_grp_i] ) ) {
      Cell_Active_Updated[i, active_grp_i] <- FALSE
      cell_final_efficacy[i, active_grp_i] <- TRUE
      Status[i, active_grp_i] <- "Close Efficacy"
    }
  }
  
  # Mark cells that remain open for recruitment.
  Status[Cell_Active_Updated & is.na(Status)] <- "Active"
  
  return(
    list(
      Cell_Active_Updated = Cell_Active_Updated,
      Cell_Final_Efficacy_Updated = cell_final_efficacy,
      Cell_Final_Futility_Updated = cell_final_futility,
      Status = Status
    )
  )
}

# Apply final futility and efficacy rules to cells that remain active.
final_decision <- function(
    theta_posterior, 
    cell_active, 
    cell_final_efficacy, 
    cell_final_futility,
    efficacy_bound, 
    futility_bound,
    efficacy_threshold, 
    futility_threshold
){
  
  Nrow<-nrow(cell_active)
  Ncol<-ncol(cell_active)
  
  Status<-matrix(" ",nrow=Nrow,ncol=Ncol)
  
  Status[!cell_active] <- "Closed" 
  cell_active_vec<-as.vector(t(cell_active)) 
  
  if(length(futility_bound)==1){
    P_futility<-apply(theta_posterior<futility_bound ,2, mean )
  } else{
    futility_bound_cell<-convert2cell(futility_bound, Nrow=Nrow, Ncol = Ncol)
    futility_bound_cell_mat<-t(replicate(nrow(theta_posterior), futility_bound_cell))
    if(!all(dim(futility_bound_cell_mat) == dim(theta_posterior))){
      message("dimension mismatch for futility bound")
    }
    P_futility<-apply(theta_posterior<futility_bound_cell_mat, 2, mean)
  }
  # Apply the final futility rule to active cells.
  Futility_Decision_vec<-rep(NA, Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Futility_Decision_vec[x] <- P_futility[x]>futility_threshold
  }
  Futility_Decision<-matrix(Futility_Decision_vec, nrow=Nrow, ncol=Ncol, byrow=TRUE)
  
  Status[Futility_Decision]<-"Conclude Futility"     
  cell_final_futility[Futility_Decision]<-TRUE
  
  # Calculate posterior efficacy probabilities.
  if(length(efficacy_bound)==1){
    P_efficacy<-apply(theta_posterior>efficacy_bound, 2, mean )
  } else{
    efficacy_bound_cell<-convert2cell(efficacy_bound,Nrow = Nrow, Ncol=Ncol)
    efficacy_bound_cell_mat<-t(replicate(nrow(theta_posterior), efficacy_bound_cell))
    
    if(!all(dim(efficacy_bound_cell_mat) == dim(theta_posterior))) {
      message("dimension mismatch for efficacy bound")
    }
    P_efficacy<-apply(theta_posterior>efficacy_bound_cell_mat, 2, mean)
  }
  
  
  Efficacy_Decision_vec <- rep(NA, Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Efficacy_Decision_vec[x] <- P_efficacy[x]>efficacy_threshold
  }
  Efficacy_Decision<-matrix(Efficacy_Decision_vec, nrow=Nrow, ncol=Ncol, byrow=TRUE)
  
  Status[Efficacy_Decision]<-"Conclude Efficacy" 
  cell_final_efficacy[Efficacy_Decision]<-TRUE
  
  
  return(
    list(
      Cell_Final_Efficacy_Updated = cell_final_efficacy,
      Cell_Final_Futility_Updated = cell_final_futility,
      Prob_Efficacy       = matrix(P_efficacy,nrow=Nrow,ncol=Ncol,byrow=TRUE),
      Prob_Futility       = matrix(P_futility,nrow=Nrow,ncol=Ncol,byrow=TRUE),
      Efficacy_Decision   = Efficacy_Decision,
      Futility_Decision   = Futility_Decision,
      Status = Status
    )
  )
}






# Declare subgroup-level success when every truly sensitive subgroup has
# at least one truly effective arm that is declared efficacious.
trial_success_subgroup <- function(
    efficacy_conclusion, 
    true_effective, 
    N_subgrp, 
    N_arm
){
  efficacy_conclusion_Mat <- matrix(efficacy_conclusion, ncol=N_subgrp, nrow=N_arm, byrow=TRUE)
  true_effective_Mat <- matrix(true_effective, ncol=N_subgrp, nrow=N_arm, byrow=TRUE)
  
  true_sensitive_subgroups <- which(colSums(true_effective_Mat) > 0)
  successful_detect <- rep(NA,N_subgrp)
  
  for (k in true_sensitive_subgroups) {
    successful_detect[k]  <- any(true_effective_Mat[,k] & efficacy_conclusion_Mat[,k])
  }
  
  return(all(successful_detect,na.rm=TRUE))
}


# Declare arm-level success when every truly effective arm has at least
# one truly sensitive subgroup that is declared efficacious.
trial_success_arm <- function(
    efficacy_conclusion, 
    true_effective, 
    N_subgrp, 
    N_arm
){
  efficacy_conclusion_Mat <- matrix(efficacy_conclusion, ncol=N_subgrp,nrow=N_arm,byrow=TRUE)
  true_effective_Mat <- matrix(true_effective,           ncol=N_subgrp,nrow=N_arm,byrow=TRUE)
  
  true_useful_arms <- which(rowSums(true_effective_Mat) > 0)
  successful_detect <- rep(NA,N_arm)
  
  for (i in true_useful_arms) {
    successful_detect[i]  <- any(true_effective_Mat[i,] & efficacy_conclusion_Mat[i,])
  }
  return(all(successful_detect,na.rm=TRUE))
}



