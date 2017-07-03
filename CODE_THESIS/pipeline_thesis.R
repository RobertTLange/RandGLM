## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

# PIPELINE FOR REPRODUCTION OF FIGURES/RESULTS IN PROJECT
####################################################

rm(list=ls())
options(scipen=999)
Sys.setenv(LANG = "en")

##### SET SEED for reproducibility ####
set.seed(101, kind = NULL, normal.kind = NULL)

## Set binary variable to false if you don't want to save plots/sim csv results
save.plots=TRUE
save.sim=TRUE

####################################################
# SETUP - Load Packages and Functions
## Install all packages that are needed
list.packs <- c("ggplot2", "latex2exp", "foreach",
                "doParallel", "corpcor", "OPDOE",
                "Matrix", "MASS", "mvtnorm", "reshape2")

new.packs <- list.packs[!(list.packs %in% installed.packages()[,"Package"])]
if(length(new.packs)) install.packages(new.packs)

## Load in libraries
sapply(list.packs, require, character.only=T)

## Load in basic functions (euclidean matrix row norm/multiplot/Hadamard Gen.)
source("A_Basic.R")
## Load in functions to compute the random sampling estimator for LS
source("B_RandLS.R")
## Load in functions to compute the random sampling estimator for GLMs
source("C_RandGLM.R")
## Load in function that aligns subplits nicely
source("multiplot.R")

####################################################
# REPRODUCTION - Step by step reproduction of figures in thesis
## Figure 1: Randomized Hadamard Transform and Leverage Score Uniformization
n = 128; d = 2
system.time(source("F1_visualization_transform.R"))
print("Finished Reproducing Fig 1")
## Figure 2: Leverage Score Sampling Distribution Approximation
n=4096; d=2; eps=0.2
system.time(source("F2_lev_score_approx.R"))
print("Finished Reproducing Fig 2")
## Figure 3: Quality of Approximation - Random Sampling LS Estimator
n=1000; d=5; iter=1000; eps=0.01
system.time(source("F3_sim_ls.R"))
print("Finished Reproducing Fig 3")
## Figure 4: Euclidean Norm Error for the Random Sampling IWLS Estimator
n=1000; d=5; iter=1000; eps=0.01; delta=1e-4; max.iter=50
system.time(source("F4_sim_glm.R"))
print("Finished Reproducing Fig 4")
## Figure 5: Trace of Estimator based on Random Sampling IWLS Estimators
n=3000; d=5; eps=0.01; delta=1e-4; max.iter=50
system.time(source("F5_glm_trace.R"))
print("Finished Reproducing Fig 5")
