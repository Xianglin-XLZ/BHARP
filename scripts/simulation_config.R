
# ---- Setup ----

N_sim      <- 500     # number of simulated trials
N_subgrp   <-10        # number of subgroups
N_arm      <-1         #number of arms


# ---- Scenario ----

scenario_list<-list(
  S_A = list(
    name = "S_A",
    Theta = matrix(c(rep(0.00,10)), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_B = list(
    name = "S_B",
    Theta = matrix(c(-.08, -.08, -.04, -.04 , 0 ,
                     0, 0.04, 0.04, 0.08, 0.08 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_C = list(
    name = "S_C",
    Theta = matrix(c(rep(0.00,7), rep(1.3,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes= matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_D = list(
    name = "S_D",
    Theta = matrix(c(-.08, -.04, -.04, 0 , 0.04 ,
                     0.04, 0.08, 1.3-0.08, 1.3-0.04 , 1.3 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(35,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_E = list(
    name = "S_E",
    Theta = matrix(c(0, 0, 0, 0, 0 ,
                     0, 0, 0.65, 0.65, 0.65 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(70,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_F = list(
    name = "S_F",
    Theta = matrix(c(-.08, -.04, -.04, 0 , 0.04 ,
                     0.04, 0.08, 0.65-0.04, 0.65 , 0.65+0.04 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(70,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_G = list(
    name = "S_G",
    Theta = matrix(c(-.08, -.04, -.04, 0 , 0.04 ,
                     0.04, 0.08, 0.65-0.04, 0.65 , 0.65+0.04 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(c(70,70,70,70,70,
                          70,20,20,70,70), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_H = list(
    name = "S_H",
    Theta = matrix(c(0, 0, 0, 0, 0 ,
                     0.65, 0.65, 0.65, 0.65, 0.65 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(70,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_I = list(
    name = "S_I",
    Theta = matrix(c(-.08, -.04, 0, 0.04 , 0.08 ,
                     0.65-.08,0.65-.04, 0.65, 0.65+0.04, 0.65+0.08 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(rep(70,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_J = list(
    name = "S_J",
    Theta = matrix(c(-.08, -.04, 0, 0.04 , 0.08 ,
                     0.65-.08,0.65-.04, 0.65, 0.65+0.04, 0.65+0.08 ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
    cell_sizes = matrix(c(70,70,70,70,20,
                          20,70,70,70,70), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_K = list( name = "S_K", 
              Theta = matrix(c(rep(0.00,4), rep(0.65,3), rep(1.30,3) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE), 
              cell_sizes = matrix(rep(70,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ),
  S_L = list( name = "S_L", 
              Theta = matrix(c(rep(0.00,7), rep(0.65,2), rep(1.30,1) ), nrow=N_arm,ncol=N_subgrp,byrow=TRUE),
              cell_sizes = matrix(rep(70,10), nrow = N_arm, ncol = N_subgrp,byrow=TRUE) ))





# ---- Hyperparameters----


hyperparam_list<-list(a_cell=5,b_cell=6,
                      a_between=4,b_between=4,
                      a_within=70,b_within=0.10^2*71)



BLAST_hyperparam<-hyperparam_list

BHM_hyperparam<-list(a_cell=5,b_cell=6,
                     a_between= 4,b_between=4)

IND_hyperparam<-list(a_cell=5,b_cell=6,
                     c_thetaind=0,p_thetaind=0.1)
