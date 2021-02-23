# -------------------------------------------------------------------------------------------------
#  Construct HMM for Model-1: 3 state model on 3 features ('hf', 'p5', 'sh')
#  Par0 are optimized by identifying starting values resulting in best fit model from 25 iterations
#  Species is included as a fixed effect on transition probabilities.

#  Script by T. Michelot and M. Conners for Conners et al 2021:
# "Hidden Markov models identify major movement modes in accelerometer and magnetometer data from four albatross species." Movement Ecology
# Contact M. Conners (connersm@gmail.com) for correspondence
# -------------------------------------------------------------------------------------------------
  

library(momentuHMM)
library(tictoc)
library(tidyverse)

# -------------------------------
# Prep dataset for HMM analysis
# -------------------------------

set.seed(3243) # for reproducible random number generation

# prep dataframes
var_names <- c("hf", "p5", "sh") # select final feature set (Model-1)
m1<-m %>% dplyr::select(all_of(var_names),"ID", "spp")
m1$ID<-as.factor(m1$ID)
m1$spp<-as.factor(m1$spp)

# convert to mementuHMM data object
data <- prepData(data = m1, coordNames = NULL)

# -------------------------------------------------
# Set model parameters - feature distributions
# -------------------------------------------------
dist <- list(hf = "weibull", p5 = "weibull", sh="weibull")
# explore shape and scale of Weibull distribution: 
# test<-rweibull(100000, shape=45, scale = 7)
# hist(test,1000, xlim=c(0,10),ylim=c(0,500))

# ----------------------------------------------------
# Run 25 iters to find optimal starting values of Par0
# ----------------------------------------------------
n_iter <- 25


tic("model1-opt: 25 iters for opitmizing Par - 3feat, 3state, spp fixed cov")
mod_list <- pbmclapply(as.list(1:n_iter), function(i) {
  
  # provide ranges (min max of Par0 for each feature for each state)
  Par0 <- list( hf     = runif(6, 
                                   c(13,.6,.85,   2.5,0.5,.9 ),  # c(state1_min_shape, state2_min_shape, state3_min_shape  state1_min_scale, state2_min_scale, state3_min_scale)
                                   c(17,.8,1.0,   3.1,0.7,1.2 )),  # c(state1_max_shape, state2_max_shape, state3_max_shape  state1_max_scale, state2_max_scale, state3_max_scale)
                p5     = runif(6,
                                   c(6.8,4.5,13, 1.2,1.6,.9),
                                   c(7.8,5.0,19, 1.6,1.9,1.1)),
                sh     = runif(6,
                                   c(1.7,2.0,1.1,  25,35,8 ),
                                   c(3.1,3.0,2.1,  30,40,11))) 
  
  # Run iter of HMM 
  mod <- fitHMM(data, 
                nbStates = 3,  # 3 states
                dist = dist,   # feature distribution definitions
                formula=~spp,  # species included as a fixed effect on transition probabilities
                Par0 = Par0)   # starting values
  return(mod)
  
  # message_parallel(paste0("iter_",i,Sys.time()))
  
}, mc.cores = 6)
toc()

# ----------------------------------------------------
# Identify best model and save
# ----------------------------------------------------
all_nllk <- unlist(lapply(mod_list, function(v) v$mod$minimum))
whichbest <- which.min(all_nllk) # Index of best fitting model (smallest negative log-likelihood)
model1_opt <- mod_list[[whichbest]] # Best fitting model
data$state <- viterbi(model4_opt) # Use Viterbi algorithm for behavioral sequence

save(mod_list, all_nllk, model1_opt, file=paste0(~,'/models_optim/model1_optim_20210119.Rdata'))
write.csv(data,file=paste0(~,'/models_optim/model1_optim_20210119.csv'),row.names=FALSE)


###### Plot state distributions and pseudo-residuals
layout(mat = matrix(c(1, 2, 3), nrow = 1))
plot(model1_opt, breaks = 30, ask = FALSE)
plotPR(model1_opt, lag.max = NULL, ncores = 1)



