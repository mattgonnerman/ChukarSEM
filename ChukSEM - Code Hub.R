### NDOW Upland Game Bird SEM
# Load Packages and Set Data Subset Info
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "LaplacesDemon", "parallel", "coda", "MCMCvis"),
       require, character.only = T)

#Last year from which data will be used
cutoff.y <- 2014 

#What is the last year of Chukar site abundance (Should be 1 + cutoff.y)
cutoff.y.chuk <- 2015 

#Last year to predict 
final.y <- year.hold <- 2015 

n.add.y <- final.y - cutoff.y
cut <- length(1976:cutoff.y) + n.add.y #Reference used to subset dataframes later

### Run Data prep 
source("./ChukSEM - Data Prep.R")

### Run Hunter Effort Solo Model to get estimates of H to create spline inputs
# source("./ChukSEM - Hunter Effort Model.R")

### Load Model Code
source("./ChukSEM - Model Only - HNC w Covs.R")

### Run Full Model to Produce Estimates
source("./ChukSEM - Nimble Prep - HNC w Covs.R")

### Save estimates and make preliminary Graphs
source("./ChukSEM - Estimate Check.R")

### Prepare estimates for shiny app
source("./ChukSEM - ShinyDataPrep.R") 
