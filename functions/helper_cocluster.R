# -------------------------------------------------------------------------
# Posterior Co-Clustering Helpers
#
# This file defines functions for calculating posterior co-clustering
# probabilities from subgroup allocation samples.
#
# Functions:
#   CoClusterProb():
#     Calculates pairwise co-clustering probabilities for one arm.
#
#   ArmListCoCluster():
#     Calculates co-clustering probabilities separately for each arm.
#
# These functions are shared by the BHARP and BLAST analyses.
# -------------------------------------------------------------------------




# to calculate co-cluster probability for one arm
CoClusterProb<-function(iz){
  s <- ncol(iz)     
  Prob_matrix <- matrix(NA, nrow = s, ncol = s)   
  for (i in seq_len(s-1L)) {
    for (j in seq.int(i + 1L, s)) {
      Prob_matrix[i, j] <- mean(iz[, i] == iz[, j]) 
    }
  }
  diag(Prob_matrix)<-1
  return(round(Prob_matrix,4))
}

# Calculate co-clustering probabilities for all arms
ArmListCoCluster<-function(combined_iz,Nsubgrp,Narm){
  
  list_iz<- list()  
  for (i in 1:Narm){
    col_start <- (i - 1) * Nsubgrp + 1  
    col_end <- col_start + Nsubgrp - 1  
    list_iz[[i]] <- combined_iz[, col_start:col_end, drop = FALSE]
  }  
  probability_list <- lapply(list_iz, CoClusterProb)
  return(probability_list)
}




