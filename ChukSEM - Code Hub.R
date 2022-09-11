### NDOW Upland Game Bird SEM
# Load Packages and Set Data Subset Info
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "LaplacesDemon", "parallel", "coda", "MCMCvis"),
       require, character.only = T)

'%notin%' <- Negate('%in%')

#Last year from which data will be used
cutoff.y <-2011

#What is the last year of Chukar site abundance (Should be 1 + cutoff.y)
cutoff.y.chuk <- 2012 

#Last year to predict 
final.y <- year.hold <- cutoff.y + 3 

n.add.y <- final.y - cutoff.y

#############################################################################################
### MULTI-SPECIES MODEL
#############################################################################################
cut <- length(1976:cutoff.y) + n.add.y 

#Reference used to subset dataframes later### Run Data prep 
source("./Code/ChukSEM - Data Prep.R")

### Load Model Code
source("./Code/ChukSEM - Model Only - HNC w Covs.R")

### Run Full Model to Produce Estimates
source("./Code/ChukSEM - Nimble Prep - HNC w Covs.R")

### Save estimates and make preliminary Graphs
source("./Code/ChukSEM - Estimate Check.R")

### Prepare estimates for shiny app
source("./Code/ChukSEM - ShinyDataPrep.R")

#############################################################################################
### CHUKAR-ONLY MODEL
#############################################################################################
list.all.y <- 1990:final.y
cut <- length(1990:cutoff.y) + n.add.y #Reference used to subset dataframes later

### Load Model Code
source("./Code/Chuk Only - Model.R")

### Run Data prep 
source("./Code/Chuk Only - Data Prep.R")

### Run Full Model to Produce Estimates
source("./Code/Chuk Only - Nimble Prep.R")

### Save estimates and make preliminary Graphs
source("./Code/Chuk Only - Estimate Check.R")

### Prepare estimates for shiny app
source("./Code/Chuk Only - ShinyDataPrep.R")


#############################################################################################
### Send to shinyapps.io
# install.packages('rsconnect')
rsconnect::setAccountInfo(name='mattgonnerman',
                          token='5920CD9D8BE3FD293A8AAF8B5676ED4D',
                          secret='ECJd+QuMwbWOHMmgxbgEngS5nK3aAGdu0DJPTTkd')
rsconnect::deployApp('./NDOWSEM/')


