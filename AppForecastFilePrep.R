lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr", "tidyr"), require, character.only = T)


#Load model outputs to get beta and sigma estimates
load(file = "./www/model_output_FullModel_predict.rdata")
beta.output <- files[[2]]
rho.output <- files[[1]]

### Hunter Effort Regression Coefficients
B.spl.hunt.df <- MCMCsummary(beta.output, 'beta.spl.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         K = sub('.*\\,', '', RowID)) %>%
  mutate(K = as.numeric(str_sub(K,1,nchar(K)-1)),
         ID = "b.SPL") %>%
  dplyr::select(ID, Species, Region, K, mean, sd) %>%
  pivot_wider(names_from = K, names_prefix = "K", values_from = c(mean, sd)) %>%
  group_split(Region)

A.hunt.df <- MCMCsummary(beta.output, 'alpha.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID))%>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "b.0") %>%
  dplyr::select(ID, Species, Region, mean, sd) 
B.econ.hunt.df <- MCMCsummary(beta.output, 'beta.econ.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.ECON") %>%
  dplyr::select(ID, Species, mean, sd) 


### Hunter Effort Covariance Matrix
rho.hunt.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 5, ncol = 5))
rho.hunt.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.hunt.est)){
  rho.hunt.list[[rho.hunt.est$Region[i]]][rho.hunt.est$Species1[i], rho.hunt.est$Species2[i]] <- rho.hunt.est$Estimate[i]
}

### Total Harvest Regression Coefficients
A.harv.df <- MCMCsummary(beta.output, 'alpha.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID))%>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "alpha.harv") %>%
  dplyr::select(ID, Species, Region, mean, sd)
  
B.ws.harv.df <- MCMCsummary(beta.output, 'beta.wintsev.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.WINT") %>%
  dplyr::select(ID, Species, mean, sd)

B.pdsi.harv.df <- MCMCsummary(beta.output, 'beta.pdsi.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.PDSI") %>%
  dplyr::select(ID, Species, mean, sd) 

B.bbs.harv.df <- MCMCsummary(beta.output, 'beta.bbs.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.BBS") %>%
  dplyr::select(ID, Species, mean, sd) %>%
  
B.spl.hunt.df <- MCMCsummary(beta.output, 'beta.spl.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         K = sub('.*\\,', '', RowID)) %>%
  mutate(K = as.numeric(str_sub(K,1,nchar(K)-1)),
         ID = "b.SPL") %>%
  dplyr::select(ID, Species, Region, K, mean, sd) %>%
  pivot_wider(names_from = K, names_prefix = "K", values_from = c(mean, sd)) %>%
  group_split(Region)


### Hunter Effort Covariance Matrix
rho.harv.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 5, ncol = 5))
rho.harv.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.harv.est)){
  rho.harv.list[[rho.harv.est$Region[i]]][rho.harv.est$Species1[i], rho.harv.est$Species2[i]] <- rho.harv.est$Estimate[i]
}




# Estimate H/N
econ <- 0
bbs <- 0
awssi <- 0
pdsi <- 0


# Splines
require(splines)
bs_bbase <- function(x, xl = min(x, na.rm = TRUE), xr = max(x, na.rm=TRUE), nseg = 5, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(bs(x, knots = knots, degree = deg, Boundary.knots = c(knots[1], knots[length(knots)])), nrow = length(x))
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]

  return(bs_matrix)
}

load("model_output_HuntEff_pred.rdata")
harvonly.H <- files[[1]]
hunter.prime   <- MCMCpstr(harvonly.H, 'H')$H #Extract hunter numbers from Model1

cut <- ncol(hunter.prime)
nseg <- 10 #Number of spline segments
BM <- array(NA, dim = c(cut,nseg+3,5,2))
Z  <- array(NA, dim = c(cut,nseg+2,5,2))
D <- diff(diag(ncol(BM[,,1,1])), diff = 1)
Q <- t(D) %*% solve(D %*% t(D))

for(i in 1:5){
  for(j in 1:2){
    BM[,,i,j] <- bs_bbase(hunter.prime[i,,j], nseg = 10)
    Z[,,i,j] <-  BM[,,i,j]%*% Q
  }
}

ZZ <- Z
ZZ[is.na(ZZ)] <- 0

# Time
time <- 1:cut

BM1 <- array(NA, dim = c(cut,nseg+3,5,2))
Z1  <- array(NA, dim = c(cut,nseg+2,5,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))

for(i in 1:5){
  for(j in 1:2){
    BM1[,,i,j] <- bs_bbase(time, nseg = 10)
    Z1[,,i,j] <-  BM1[,,i,j]%*% Q1
  }
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0








# y.start.fore <- 2015
# y.stop.fore <- 2022
# 
# ###MAKE SURE THE BELOW MATCHES THE ORIGINAL MODEL RUN
# ### Run Initial Data Management
# lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
# cutoff.y <- 2016 #Last year from which data will be used
# final.y <- 2017 #Last year to predict
# 
# drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models
# 
# n.add.y <- final.y - cutoff.y
# source("./ChukSEM - Data Prep - Predict.R")
# 
# #BBS Index
# bbs.fore <- read.csv("./Data/bbs_indices.csv") %>%
#   mutate(raven = scale(raven)[,1],
#          rthawk = scale(rthawk)[,1],
#          nharrier = scale(nharrier)[,1],
#          pfalcon = scale(pfalcon)[,1]
#   ) %>%
#   filter(Year >= y.start.fore & Year <= y.stop.fore)
# 
# predict.fore1 <- rbind(bbs.fore %>% dplyr::select(year = Year, raven, nharrier),
#                       data.frame(year = (1+max(bbs.fore$Year)):y.stop.fore,
#                                  raven = NA, nharrier = NA))
# 
# #Drought Index
# eastern.pdsi.fore <- read.csv('./Data/Eastern_PDSIZ.csv') %>%
#   mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
#          Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
#   select(Year, Summer, Winter) %>%
#   filter(Year %in% y.start.fore:y.stop.fore) %>%
#   mutate(Region = "Eastern")
# 
# western.pdsi.fore <- read.csv('./Data/Western_PDSIZ.csv') %>%
#   mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
#          Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
#   select(Year, Summer, Winter) %>%
#   filter(Year %in% y.start.fore:y.stop.fore) %>%
#   mutate(Region = "Western")
# 
# pdsi.fore <- rbind(eastern.pdsi.fore, western.pdsi.fore) %>%
#   dplyr::rename(year = Year, pdsi_sum = Summer, pdsi_win = Winter, region = Region)
# 
# predict.fore2 <- merge(predict.fore1, pdsi.fore, by = "year", all.x = T) %>%
#   arrange(year, region)
# 
# #Unemployment
# une.mean <- attr(scale(unemployment$Rate), "scaled:center")
# une.sd <- attr(scale(unemployment$Rate), "scaled:scale")
# unemployment.fore <- read.csv("./Data/NEvadaUnemploymentUSBL.csv", colClasses = c("integer", "character", rep("numeric", 2), rep("integer", 3), "numeric")) %>%
#   filter(Period == "Sep") %>%
#   filter(Year %in% y.start.fore:y.stop.fore) %>%
#   mutate(Rate = unemployment/labor.force) %>%
#   mutate(Rate.Z = ifelse(is.na(Rate), NA, (Rate - une.mean)/une.sd)) %>%
#   dplyr::select(year = Year, une = Rate.Z)
# predict.fore3 <- merge(predict.fore2, unemployment.fore, by = "year", all.x = T)
# 
# 
# #Number of Resident Licenses Sold
# res.mean <- attr(scale(general_nv$Resident.Licenses), "scaled:center")
# res.sd <- attr(scale(general_nv$Resident.Licenses), "scaled:scale")
# general_nv.fore <- subset(general, ST == 'NV' & Year %in% y.start.fore:y.stop.fore) %>%
#   mutate(Resident.Licenses.Z = ifelse(is.na(Resident.Licenses), NA, (Resident.Licenses - res.mean)/res.sd)) %>%
#   dplyr::select(year = Year, res = Resident.Licenses.Z)
# predict.fore4 <- merge(predict.fore3, general_nv.fore, by = "year", all.x = T)
# 
# #
# gas.march.fore <- read.csv('./Data/NV_Gas.csv') %>%
#   mutate(Date = as.Date(M.Y, format = '%m/%d/%Y')) %>%
#   mutate(Year = lubridate::year(Date),
#          Month = lubridate::month(Date)) %>%
#   filter(Month == 5) %>%
#   arrange(Year) %>%
#   filter(Year %in% y.start.fore:y.stop.fore) %>%
#   dplyr::select(Year, Gas.May = Dollars.per.Gallon)
# 
# econ_data.fore <- read.csv('./Data/economic_data.csv') %>%
#   dplyr::rename(PDI = Per.Capita.Personal.Disposable.Income, Gas.Unk = Gas.Price..July.Oct.) %>%
#   merge(., gas.march.fore, by = "Year", all = T) %>%
#   filter(Year %in% y.start.fore:y.stop.fore) %>%
#   mutate(PDI = PDI/10000)
# 
# predict.fore$PDI <- econ_data$PDI
# predict.fore$GAS <- econ_data$Gas.May
