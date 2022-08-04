### Set Seed For Model Construction
set.seed(2)

### Load Packages
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "coda", "MCMCvis", "parallel", "coda"),
       require, character.only = T)


#################
### Data Prep ###
#################
### Hunter Numbers for each Species*Year
#Load all species harvest data
harvest_data <- read.csv('./Data/all_species_harvest_data.csv') %>%
  mutate(Year = as.numeric(as.character(Year)),
         Animals = as.numeric(as.character(Animals)),
         Hunters = as.numeric(as.character(Hunters)),
         H.Days = as.numeric(as.character(H.Days))) %>%

#Add species codes
harvest_data$Species_Code <- recode(harvest_data$Species, 
                                    `SAGE_GROUSE` = "SAGR",
                                    `SOOTY_GROUSE` = "BLGR",
                                    `RUFFED_GROUSE` = "RUGR",
                                    `CHUKAR` = "CHUK",
                                    `HUNGARIAN_PART` = 'HUPA',
                                    `CALIFORNIA_QUAIL` = 'CAQU',
                                    `MOUNTAIN_QUAIL` = 'MOQU',
                                    `PHEASANTS` = 'PHEA',
                                    `RABBIT` = 'RABB',
                                    `DOVE` = 'DOVE',
                                    `DUSKY_GROUSE` = 'BLGR',
                                    `GAMBEL'S_QUAIL` = 'GAQU'
)

#Subset data from Eastern section
hunters_wide <- cast(harvest_datae, Species_Code  ~ Year, value = 'Hunters', fun.aggregate = 'sum', na.rm = TRUE)
hunters_wide[hunters_wide==0] <-NA #Change 0 values to NA
species <- hunters_wide$Species_Code #Extract vector of species codes
hunters_wide$Species_Code  <- NULL #Remove species code column from df
hunters_wide<- data.matrix(hunters_wide) #convert df to matrix

#Same thing as previous but for western section
harvest_dataw <- subset(harvest_data, grepl('Western', Section))
hunters_widew <- cast(harvest_dataw, Species_Code  ~ Year, value = 'Hunters', fun.aggregate = 'sum', na.rm = TRUE)
hunters_widew[hunters_widew==0] <-NA
speciesw <- hunters_widew$Species_Code 
hunters_widew$Species_Code  <- NULL
hunters_widew<- data.matrix(hunters_widew)

hunters_widew <- hunters_widew[-6, ]#Remove MOQU

#Combine datasets from East and West sections into arrays
hunters <- abind(hunters_widew, hunters_wide, along = 3)

#Remove additional rows associated with ?
hunters <-hunters[-c(4,8),,]

#Change closed season values to NA
hunters[7,10,] <- NA # Season closed 

### Palmer Drought Severity Index (PDSI)
#Load year/month specific Eastern PDSI values (Z standardized?)
eastern_pdsi <- read.csv('./Data/Eastern_PDSIZ.csv') %>%
  filter(Year > 1974) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (lag(Nov) + lag(Dec) + Jan + Feb + Mar)/5) %>%
  select(Year, Summer, Winter) %>%
  mutate(Region = "Eastern")

western_pdsi <- read.csv('./Data/Western_PDSIZ.csv') %>%
  filter(Year > 1974) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (lag(Nov) + lag(Dec) + Jan + Feb + Mar)/5) %>%
  select(Year, Summer, Winter) %>%
  mutate(Region = "Western")

pdsi_df <- rbind(eastern_pdsi, western_pdsi) %>%
  pivot_wider(names_from = Region, values_from = c(Winter, Summer))

#Average breeding and winter DPSI values for each section
pdsi <- as.matrix(pdsi_df[-1,] %>% select(Summer_Eastern, Summer_Western)) 
wpdsi <- as.matrix(pdsi_df[-1,] %>% select(Winter_Eastern, Winter_Western))

cut <- 42 #Number of years to keep (1976-2017)

### Economic Metrics ###
#Unemployment data, 1976 - 2021 (Bureau of Labor, https://www.bls.gov/regions/west/nevada.htm#eag)
unemployment <- read.csv("./Data/NEvadaUnemploymentUSBL.csv", colClasses = c("integer", "character", rep("numeric", 2), rep("integer", 3), "numeric")) %>%
  filter(Period == "Sep") %>%
  filter(Year > 1975) %>%
  mutate(Rate = unemployment/labor.force)
une <- scale(unemployment$Rate)

# unemployment <- c(8.825,6.766666667,4.416666667,4.866666667,6.316666667,7.316666667,10.05,9.841666667,7.658333333,7.475,6.383333333,6.2,5.241666667,4.658333333,4.7,5.9,6.8,6.9,6.241666667,5.566666667,5.058333333,4.375,4.175,3.966666667,4.116666667,5.191666667,5.683333333,5.283333333,4.458333333,4.125,4.066666667,4.591666667,6.883333333,11.725,13.73333333,13.31666667,11.60833333,9.966666667,8.158333333,6.85,5.808333333,5.016666667,4.408333333,3.9,13.01666667,8.333333333)
# unemployment_cut <- unemployment[1:(44+1)]
# une <- scale(unemployment_cut)

#Number of Resident Licenses Sold
general <- read.csv(file = './Data/overall_hunting_data.csv')
general_nv <- subset(general, ST == 'NV' & Year > 1975)
res <- scale(general_nv$Resident.Licenses)[,1]

#Gas Prices
econ_data <- read.csv('./Data/economic_data.csv')
econ_data_sub <- subset(econ_data, Year > 1975)
PDI <- (econ_data_sub$Per.Capita.Personal.Disposable.Income/10000)
GAS <- econ_data_sub$Gas.Price..July.Oct.

#Objects not used but values are represented in model
mean_ratio <- 2.581635
sd_ratio <- 0.8894599

#Calculate ratio of disposable income to gas price
ratio <- PDI/GAS

#Accumulated Winter Season Severity Index
#https://mrcc.purdue.edu/research/awssi/indexAwssi.jsp#
awssi.df <- read.csv("./Data/Nevada AWSSI.csv") %>%
  mutate(Start = as.Date(Start, "%m/%d/%Y")) %>%
  mutate(Year = lubridate::year(Start)) %>%
  filter(Region != "Southern", Year > 1975) %>%
  group_by(Region, Year) %>%
  summarise(AWSSI = mean(AWSSI)) %>%
  pivot_wider(names_from = Region, values_from = AWSSI)

awssi <- t(awssi.df[2:3])


########################
### Prelim Model Run ###
########################
### This initial model run provides information for the spline in a second run ###

### Prepare Data Inputs
## Splines
require(splines)
### Function that Constructs B-Spline Base
### from? = https://github.com/andrewcparnell/jags_examples/blob/master/R%20Code/jags_spline.R
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

nseg <- 10
BM1 <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z1  <- array(NA, dim = c(cut+4,nseg+2,7,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))
time <- 1:(cut+4)
for(i in 1:7){
  for(j in 1:2){
    BM1[,,i,j] <- bs_bbase(time, nseg = 10)
    Z1[,,i,j] <-  BM1[,,i,j]%*% Q1
    
  }
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0

### Set Initial Values
nu = n.species + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.species, n.species)
for (j in 1:n.species){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)

Hi <- hunters + 2500
Hi[,-1,] <- NA
Hi[7,10,] <- rpois(2,colMeans(hunters[7,-10,]))
Hi <- abind(Hi[,1:cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi <- abind(hunters,hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi[7,10,] <- Hi[7,10,]  

rho.hunt.init <- diag(7)
PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), NA)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, NA)

### Define Constants
require(LaplacesDemon)
n.species<- dim(hunters)[1]
sig = rgamma(n.species,1,1)
Lambda = diag(sig)
I = diag(n.species) #identity matrix

# Package as List for NIMBLE
data <- list(une = une[,1], #unemployment
             PDI = PDI, #personal disposable income
             GAS = GAS, #gas prices
             n.hunt = hunters, #number of hunters (array = region)
             Z.hunt = ZZ1, #for spline
             res= scale(res)[,1], #residential license sales
             wpdsi = data.matrix(abind(wpdsi,matrix(NA,1,2), along = 1))) #drought metric

constants <- list(n.species = 7, #number of species
                  K= 12, #number of knots for spline
                  n.region = 2, #East versus West
                  n.year = 2017-1975, #Number of Years
                  mu.hunt = rep(0, 7), #mean change in harvest (0 centered)
                  era = c(rep(1,19),rep(2, 27)), #Groupings for change in gas prices 
                  mean.H = apply(hunters, c(1,3), mean, na.rm = TRUE), # Mean Harvest
                  sd.H = apply(hunters, c(1,3), sd, na.rm = TRUE), # SD Harvest
                  I = abind(I,I,along = 3)) #Identity Matrix

initsFunction <- function() list(  #Dan's
  GAS = GAS.inits,
  PDI = PDI.inits,
  sig.pdi = 1, 
  sig.gas = 1,
  beta.drought2 = matrix(0, 7, 2), 
  mu.drought2 = c(0,0), 
  sig.drought2 = c(1,1),
  beta.income = matrix(0, 7, 2), 
  mu.incom = c(0,0), 
  sig.incom= c(1,1),
  beta.jobs = matrix(0, 7, 2), 
  mu.jobs = c(0,0), 
  sig.jobs = c(1,1), 
  sig.wpdsi = c(1,1), 
  sig.une = 1,
  
  b0.pdi = 2.9, 
  bt.pdi = 0, 
  b0.gas = 0.9,  
  ar1 = 0, 
  bt.gas = c(0,0),
  Q = abind(Q,Q,along = 3),
  Sigma = abind(Sigma,Sigma,along = 3),
  P = abind(P,P,along = 3),
  Lambda = abind(Lambda,Lambda,along = 3),
  Delta = abind(Delta,Delta,along = 3),
  rho = abind(diag(n.species),diag(n.species),along = 3),
  
  H = Hi, 
  x = zi,
  sig.trend = matrix(1, ncol = 2, nrow = 7),
  beta.general = matrix(0, ncol = 2, nrow = 7), 
  mu.gen = c(0,0), 
  sig.gen = c(1,1),
  lbo1 = matrix(0, ncol = 2, nrow = 7),
  
  lambda1 = array(1, dim = c(7,cut+3,2) ),
  log.r1  = array(0, dim = c(7,cut+3,2) ),
  hunt.eps = array(rnorm(7*2*cut+4,0,0.1),dim = c(7,2,cut+4)),
  
  wpdsi = matrix(0,46,2),
  WPDSI = matrix(0,46,2),
  une = rep(0,46),
  UNE = rep(0,46)
)


inits <- initsFunction()

### The Model
code <- nimbleCode( {
  #Personal Disposable Income
  b0.pdi ~ dnorm(0, 0.001) #intercept
  bt.pdi ~ dnorm(0, 0.01) #year
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1 ~ dunif(-1,1) #Autoregressive parameter
  
  pdi.trend[1] <- b0.pdi + bt.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    mu.pdi[t] <- pdi.trend[t] + ar1 * ( PDI[t-1] - pdi.trend[t-1] )
    pdi.trend[t] <- b0.pdi + bt.pdi * t
  } #t
  
  #Gas Prices
  b0.gas ~ dnorm(1.5, 1)
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  bt.gas[1] ~ dnorm(0, 0.01)
  bt.gas[2] ~ dnorm(0, 0.01)
  
  #Relative Cost of Gas
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(b0.gas + bt.gas[era[t]]*t, sd = sig.gas),0,)
    REL.COST[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  } #t
  
  #Unemployment Rate
  sig.une ~ dunif(0,5)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
  } #t
  
  #Drought Index
  for(r in 1:n.region){
    sig.wpdsi[r] ~ dunif(0,5)
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
    } #t
  } #r
  
  
  for(r in 1:n.region){
    mu.income[r] ~ dnorm(0, 0.01)
    sig.income[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.license[r] ~ dnorm(0, 0.01)
    sig.license[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.drought[r] ~ dnorm(0, 0.01)
    sig.drought[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.jobs[r] ~ dnorm(0, 0.01)
    sig.jobs[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    
    ### Coefficient Values
    for(s in 1:n.species){
      sig.spline.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      alpha.hunter[s,r] ~ dnorm(5, sd = 1)
      for(k in 1:K){
        beta.spline.hunt[s,r,k] ~ dnorm(0, sd = sig.spline.hunt[s,r])
      } #k
      beta.drought[s,r] ~ dnorm(mu.drought[r], sd = sig.drought[r])
      beta.jobs[s,r] ~ dnorm(mu.jobs[r], sd = sig.jobs[r])
      beta.license[s,r] ~ dnorm(mu.license[r], sd  = sig.license[r])
      beta.income[s,r] ~ dnorm(mu.incom[r], sd  = sig.incom[r])
      
      ### Regression for Number Hunters
      for(t in 1:n.year){
        pred1[s,r,t] <- inprod(beta.spline.hunt[s,r,1:K], Z.hunt[t,1:K,s,r])
        
        mu.hunter[s,r,t] <- alpha.hunter[s,r] + #Intercept
          inprod(beta.spline.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) + #Spline smoothing terms
          beta.license[s,r] * res[t] + #Hunting Licences Sold
          beta.income[s,r] * rel.cost[t] + #PDI/Gas Price
          beta.drought[s,r] * wpdsi[t,r] + #Drought
          beta.jobs[s,r] * une[t] #Unemployment
      } #t
      
      #Quantify annual change in number of hunters
      for(t in 2:n.year){
        lambda.hunter[s,t-1,r] <- H[s,t,r]/H[s,t-1,r]
      }#t
    } #s
    
    ### Define Variance-Covariance Matrix for Number of Hunters
    Q.hunt[1:n.species,1:n.species,r] ~ dinvwish(S = I.hunt[1:n.species,1:n.species,r], df = n.species + 1) 
    
    #Variance
    for(s in 1:n.species){
      sig.hunt[s,r] ~ dgamma(1,1)
      Lambda.hunt[s,s,r] <- sig.hunt[s,r]
      Delta.hunt[s,s,r] <- pow(Q.hunt[s,s,r], -0.5)
    } #s
    #Covariance
    for (s1 in 2:n.species){
      for (s2 in 1:(s1-1)){
        Lambda.hunt[s1,s2,r] <- 0
        Delta.hunt[s1,s2,r] <- 0
      } #j
    } #i
    
    P.hunt[1:n.species,1:n.species,r] <- Delta.hunt[1:n.species,1:n.species,r] %*% Q.hunt[1:n.species,1:n.species,r] %*% Delta.hunt[1:n.species,1:n.species,r]
    
    #Covariance matrix for multivariate normal dist. describing number of hunters by species
    Sigma.hunt[1:n.species,1:n.species,r] <- Lambda.hunt[1:n.species,1:n.species,r] %*% P.hunt[1:n.species,1:n.species,r] %*% Lambda.hunt[1:n.species,1:n.species,r]
    
    #Calculate Correlation Coefficient
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        #Covariance/(sqrt(SD1)*sqrt(SD2))?
        rho[s1,s2,r] <- Sigma[s1,s2,r]/sqrt(Sigma[s1,s1,r] * Sigma[s2,s2,r])
      } #s2
    } #s1
    
    ### Process determining observed number of hunters
    for(t in 1:n.year){
      hunt.eps[1:n.species,r,t] ~ dmnorm(mu.hunter[1:n.species,r,t], cov = Sigma.hunt[1:n.species,1:n.species,r])
    } #t
    
    # Number of Hunters
    for(s in 1:n.species){
      for(t in 1:(n.year)){
        log(H[s,t,r]) <- hunt.eps[s,r,t] #Log link to restrict H above 0
        n.hunt[s,t,r] ~  dpois(H[s,t,r]) #Observed number of hunters species*year
      } #t
    } #s
  }#r
})

### Set Parameter Monitors
pars1 <- c('b0')

### Model Test
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c('sig',
                      'pred1',
                      'mu',
                      'REL.COST',
                      'mu.pdi',
                      'une',
                      'wpdsi',
                      'beta.trend'))
model_test$initializeInfo()
model_test$calculate()