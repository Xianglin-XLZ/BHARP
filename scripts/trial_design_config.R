
N_sim      <- 500      # number of simulated trials
N_subgrp   <- 6        # number of subgroup; k
N_arm      <- 3        # number of arm; i
L          <- 4        # number of analysis; l

# Probability thresholds
P_fut <- rep(0.5,L)      # Futility thresholds
P_eff <- c(0.95,0.95,0.95,0.95)     # Efficacy thresholds
# boundaries
x_fut <- 0.15                     # futility   boundaries; MCID 
x_eff <- 0.5            # efficacy boundaries

ThetaVec<-c(.1 , .1 , .1, .1 , .1 , .1,
            .1, .4 , .7,  .13,  .43, .73,
            .6, .65 , .7 ,  1.2, 1.22, 1.25)


Theta<-matrix(ThetaVec,nrow = N_arm,ncol=N_subgrp,byrow = T)    #scenario; true cell mean

TotalSampleSize<-c(1200,1500,1800,2100)    #total sample size at each analysis time point

hyperparam_list<-BLAST_hyperparam<-list(a_cell=5,b_cell=6,
                      a_between=4,b_between=0.8*5,
                      a_within=150,b_within=0.10^2*151)

BLAST_hyperparam<-hyperparam_list

IND_hyperparam<-list(a_cell=5,b_cell=6,
                     c_thetaind=0,p_thetaind=0.1)

BHM_hyperparam<-list(a_cell=5,b_cell=6,
                     a_between= 4,b_between=4)