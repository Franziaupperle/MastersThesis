library(lars) 
library(scalreg)
library(hdm) 
library(Matrix)
library(glmnet) 
library(lattice)
library(ggplot2)
library(hdi)
library(caret)
library(Metrics)
library(parallel)
library(iterators)
library(foreach)
library(doParallel)
library(Matrix)
library(MASS)
library(plotrix)


### IN THIS FILE, PARAMETER SETS AS WELL AS INTIALISING CALCULATION OF METHODS

source("Methods.R") 
source("generate-sigma.R")
source("generate-parameter.R")
source("generate-data.R")
source("Oracle.R")  
source("Double.R")
source("R-Split.R")
source("PODS.R")
source("PODS-Split.R")
source("DB-stab.R")


# Choose the parameter setting

no_runs_mc = 100
B = 1000
n = 100
p = 500
beta = "approximately sparse" # "sparse", "approximately sparse", "moderately sparse" or "dense"
gamma = "dense" # "sparse" or "dense"
type = "equalcorr9" # "ind", "toeplitz9", "equalcorr9" or"equalcorr3"
Ry2 = 0.8
Rd2 = 0.8
dfmin = 10
tresh = 0.5


# Conduct Monte Carlo Simulations for each method, estimating the treatment effect and calculating performance measures

# Oracle Simulation
oracle3 = oracle.mcfun(no_runs_mc = no_runs_mc)

# Double Simulation
double6 = double.mcfun(no_runs_mc = no_runs_mc)

# PODS Simulation
pods6 = pods.mcfun(no_runs_mc = no_runs_mc)

# R-SPLIT Simulation
rsplit6 = r.split_mcfun(no_runs_mc = no_runs_mc, B = B, dfmin = dfmin)

# PODS-SPLIT Simulation
pods.spli6 = pods.split_mcfun(no_runs_mc = no_runs_mc, B = B, dfmin = dfmin)

# DB-STABILITY Simulation
db.stab6 = db.stab_mcfun(no_runs_mc = no_runs_mc, B = B, tresh = tresh)
 
