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


### IN THIS FILE, PARAMETER SETTINGS ARE CHOSEN AND SIMULATIONS ARE STARTED

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
beta = "sparse" # "sparse", "approximately sparse", "moderately sparse" or "dense"
gamma = "sparse" # "sparse" or "dense"
type = "toeplitz9" # "ind", "toeplitz9", "equalcorr9" or"equalcorr3"
Ry2 = 0.8
Rd2 = 0.5
dfmin = 5
tresh = 0.7


# Conduct Monte Carlo Simulations for each method, estimating the treatment effect and calculating performance measures

# Simulation with the Oracle 
oracle1 = oracle.mcfun(no_runs_mc = no_runs_mc)

# Simulation with Double-Selection
double1 = double.mcfun(no_runs_mc = no_runs_mc)

# Simulation with PODS
pods1 = pods.mcfun(no_runs_mc = no_runs_mc)

# Simulation with R-SPLIT 
rsplit1 = r.split_mcfun(no_runs_mc = no_runs_mc, B = B, dfmin = dfmin)

# Simulation with PODS-SPLIT 
pods.spli1 = pods.split_mcfun(no_runs_mc = no_runs_mc, B = B, dfmin = dfmin)

# Simulation with DB-STABILITY 
db.stab1 = db.stab_mcfun(no_runs_mc = no_runs_mc, B = B, tresh = tresh)
 
