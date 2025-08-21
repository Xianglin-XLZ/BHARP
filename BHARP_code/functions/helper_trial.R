# File: helper_trial.R
# Purpose: Provide helper functions to support Bayesian adaptive trial simulations, 
#          including patient recruitment and interim/final decision-making logic
# Contents:
#   - random_balanced_partition: split an integer into approximately equal parts
#   - random_proportional_partition: split an integer according to a ratio
#   - convert2cell: convert arm-/group-level thresholds to cell-level
#   - recruit_patients: simulate new patient recruitment based on allocation strategy
#   - deactivate_decision: apply interim rules to update cell activity status
#   - final_decision: apply final rules to conclude efficacy or futility
# Dependency:
#   - Requires posterior samples of treatment effects (theta_posterior)
#   - Used by trial simulation routines for patient generation and decision logic




# helper function to split an integer into balanced integers.
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
  
  # Distribute the remainder randomly
  if (remainder > 0) {
    indices <- sample(1:parts, remainder)
    partition[indices] <- partition[indices] + 1
  }
  return(partition)
}

# helper function to split an integer by allocation ratio
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



#helper function to convert grp-wise or arm-wise bounds to cell wise,by row
#if Nrow=Ncol should input complete length(number)==Nrow*Ncol
convert2cell<-function(number,Nrow,Ncol){
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

# for simulation: 
# generate dataset , recruiting up to targetN
recruit_patients <- function(dataset, cell_active, targetN,
                             n_subgrp, n_arm, Theta, 
                             method = c("balanced", "proportion"),
                             ratio = NULL) {
  
  #   - dataset: current dataset, it can be empty 
  #   - cell_active:  matrix/logical,dimension n_arm * n_subgrp, showing which cells are recruiting
  #   - targetN: target total sample size
  #   - n_arm, n_subgrp: number of arm and subgroup
  #   - method: c("balanced", "proportion") how to recruit into active cells
  #   - ratio:  if method=="proportion" , should provide a vector with same length active cells 
  
  method <- match.arg(method)
  
  # 1) calculate current sample size
  currentN <- nrow(dataset)
  if (is.null(currentN)) currentN <- 0  
  
  # 2) calculate  sample size to recruit in this stage
  newN <- targetN - currentN
  if (newN <= 0) {
    message("No new participants needed. Current N >= targetN.")
    return(dataset)
  }
  
  # 3) Find active (i,k) cells 
  #    cell_active[i,k] == TRUE means the cell is accepting new patients 
  if (!all(dim(cell_active) == c(n_arm, n_subgrp))) {
    stop("cell_active should have dim: n_arm , n_subgrp!")
  }
  # which(..., arr.ind=TRUE) a matrix of integer pairs (i,k), each row is a pair
  active_cells <- which(cell_active, arr.ind = TRUE) 
  
  n_active <- nrow(active_cells)  # number of active cells
  
  # 4) calculate new patients in each cell
  if (method == "balanced") {
    
    new_partitions <- random_balanced_partition(newN, parts = n_active)
  } else if (method == "proportion") {
    if (is.null(ratio)) {
      # if method==proportion  ratio=NULL, ratio=1,...,1
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
  
  # 5) For each cell generate a data frame: (i,k) repeat new_partitions[idx] times
  new_recruits_list <- lapply(seq_len(n_active), function(idx) {
    i_val <- active_cells[idx, 1]
    k_val <- active_cells[idx, 2]
    count <- new_partitions[idx]
    
    if (count > 0) {
      data.frame(i = rep(i_val, count),
                 k = rep(k_val, count),
                 Y=rnorm(count, Theta[i_val,k_val], sd=1))
    } else {
      # if any cell is not recruiting return an empty data frame
      data.frame(i = integer(0), k = integer(0), Y=numeric(0))
    }
  })
  
  new_recruits <- do.call(rbind, new_recruits_list)
  
  
  
  updated_dataset <- rbind(dataset, new_recruits)
  
  return(updated_dataset)
}












# efficacy bound and futility bound: defining efficacy and futility
# probability thresholds: confidence level to trigger decisions,
# close for futility and efficacy

deactivate_decision<-function(theta_posterior, cell_active, cell_final_efficacy, cell_final_futility,
                              efficacy_bound, futility_bound,
                              efficacy_threshold, futility_threshold){
  
  Nrow<-nrow(cell_active); Ncol<-ncol(cell_active)
  
  Status<-matrix(NA,nrow=Nrow,ncol=Ncol)
  Status[!cell_active]<-"Closed" #record cells that were closed in previous analysis
  
  cell_active_vec<-as.vector(t(cell_active)) #make a copy of active status as vector
  
  if(length(futility_bound)==1){
    P_futility<-apply(theta_posterior<futility_bound ,2, mean )
  } else{
    futility_bound_cell<-convert2cell(futility_bound,Nrow=Nrow,Ncol = Ncol)
    futility_bound_cell_mat<-t(replicate(nrow(theta_posterior), futility_bound_cell))
    if(!all(dim(futility_bound_cell_mat) == dim(theta_posterior))){message("dimension mismatch for futility bound")}
    P_futility<-apply(theta_posterior<futility_bound_cell_mat ,2, mean )
  }
  
  Futility_Decision_vec<-rep(NA,Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Futility_Decision_vec[x] <- P_futility[x]>futility_threshold
    if(isTRUE(Futility_Decision_vec[x])){cell_active_vec[x]<-FALSE } #close cell for futility
  }
  Futility_Decision<-matrix(Futility_Decision_vec,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  Status[Futility_Decision]<-"Close Futility"     #Record close for futility
  cell_final_futility[Futility_Decision]<-TRUE
  
  
  if(length(efficacy_bound)==1){
    P_efficacy<-apply(theta_posterior>efficacy_bound ,2, mean )
  } else{
    efficacy_bound_cell<-convert2cell(efficacy_bound,Nrow = Nrow, Ncol=Ncol)
    efficacy_bound_cell_mat<-t(replicate(nrow(theta_posterior), efficacy_bound_cell))
    if(!all(dim(efficacy_bound_cell_mat) == dim(theta_posterior))) {message("dimension mismatch for efficacy bound")}
    P_efficacy<-apply(theta_posterior>efficacy_bound_cell_mat ,2, mean )
  }
  
  
  Efficacy_Decision_vec<-rep(NA,Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Efficacy_Decision_vec[x] <- P_efficacy[x]>efficacy_threshold
    #if(isTRUE(Efficacy_Decision[x])){cell_active_vec[x]<-FALSE}  # close cell for efficacy
  }
  Efficacy_Decision<-matrix(Efficacy_Decision_vec,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  Status[Efficacy_Decision]<-"Active Efficacy" #Have evidence not close
  
  
  Cell_Active_Updated<- matrix(cell_active_vec,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  #Prob_Futility<-matrix(P_futility,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  #Prob_Efficacy<-matrix(P_efficacy,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  
  
  #make arm-wise decisions. if all remaining active cells are effective, close this arm
  for (i in which(rowSums(Cell_Active_Updated) > 0) ){
    active_grp_i<-which(Cell_Active_Updated[i,]) 
    if (all( Efficacy_Decision[i, active_grp_i] ) ) {
      Cell_Active_Updated[i, active_grp_i] <- FALSE
      cell_final_efficacy[i, active_grp_i] <- TRUE
      Status[i, active_grp_i] <- "Close Efficacy"
    }
  }
  
  #other recruiting cells without conclusion:Active
  Status[Cell_Active_Updated & is.na(Status)]<-"Active"
  
  return(list(
    Cell_Active_Updated = Cell_Active_Updated,
    Cell_Final_Efficacy_Updated = cell_final_efficacy,
    Cell_Final_Futility_Updated = cell_final_futility,
    #Prob_Efficacy       = Prob_Efficacy,
    #Prob_Futility       = Prob_Futility,
    #Efficacy_Decision   = Efficacy_Decision,
    #Futility_Decision   = Futility_Decision,
    Status = Status
  ))
}

# efficacy bound and futility bound: defining efficacy and futility
# probability thresholds: confidence level to trigger decisions,
# conclude futility and efficacy 

final_decision<-function(theta_posterior, cell_active, cell_final_efficacy, cell_final_futility,
                         efficacy_bound, futility_bound,
                         efficacy_threshold, futility_threshold){
  
  Nrow<-nrow(cell_active); Ncol<-ncol(cell_active)
  
  Status<-matrix(" ",nrow=Nrow,ncol=Ncol)
  Status[!cell_active]<-"Closed" #record cells that were closed in previous analysis
  
  cell_active_vec<-as.vector(t(cell_active)) #make a copy of active status as vector
  
  if(length(futility_bound)==1){
    P_futility<-apply(theta_posterior<futility_bound ,2, mean )
  } else{
    futility_bound_cell<-convert2cell(futility_bound,Nrow=Nrow,Ncol = Ncol)
    futility_bound_cell_mat<-t(replicate(nrow(theta_posterior), futility_bound_cell))
    if(!all(dim(futility_bound_cell_mat) == dim(theta_posterior))){message("dimension mismatch for futility bound")}
    P_futility<-apply(theta_posterior<futility_bound_cell_mat ,2, mean )
  }
  
  Futility_Decision_vec<-rep(NA,Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Futility_Decision_vec[x] <- P_futility[x]>futility_threshold
  }
  Futility_Decision<-matrix(Futility_Decision_vec,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  Status[Futility_Decision]<-"Conclude Futility"     
  cell_final_futility[Futility_Decision]<-TRUE
  
  
  
  if(length(efficacy_bound)==1){
    P_efficacy<-apply(theta_posterior>efficacy_bound ,2, mean )
  } else{
    efficacy_bound_cell<-convert2cell(efficacy_bound,Nrow = Nrow, Ncol=Ncol)
    efficacy_bound_cell_mat<-t(replicate(nrow(theta_posterior), efficacy_bound_cell))
    if(!all(dim(efficacy_bound_cell_mat) == dim(theta_posterior))) {message("dimension mismatch for efficacy bound")}
    P_efficacy<-apply(theta_posterior>efficacy_bound_cell_mat ,2, mean )
  }
  
  
  Efficacy_Decision_vec<-rep(NA,Nrow*Ncol)
  for (x in which(cell_active_vec) ){
    Efficacy_Decision_vec[x] <- P_efficacy[x]>efficacy_threshold
  }
  Efficacy_Decision<-matrix(Efficacy_Decision_vec,nrow=Nrow,ncol=Ncol,byrow=TRUE)
  Status[Efficacy_Decision]<-"Conclude Efficacy" 
  cell_final_efficacy[Efficacy_Decision]<-TRUE
  
  
  return(list(
    #Cell_Active_Updated = Cell_Active_Updated,
    Cell_Final_Efficacy_Updated = cell_final_efficacy,
    Cell_Final_Futility_Updated = cell_final_futility,
    Prob_Efficacy       = matrix(P_efficacy,nrow=Nrow,ncol=Ncol,byrow=TRUE),
    Prob_Futility       = matrix(P_futility,nrow=Nrow,ncol=Ncol,byrow=TRUE),
    Efficacy_Decision   = Efficacy_Decision,
    Futility_Decision   = Futility_Decision,
    Status = Status
  ))
}




