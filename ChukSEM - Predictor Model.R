### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2016 #Only need to change this to adjust the number of years

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- 2017-cutoff.y
source("./ChukarSEM - 1 Data Prep - Predict.R")

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Predictors
  #Personal Disposable Income
  alpha.pdi ~ dnorm(0, sd = 100) #intercept
  beta.t.pdi ~ dnorm(0, sd = 100) #year
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1.pdi ~ dunif(-1,1) #Autoregressive parameter
  
  pdi.trend[1] <- alpha.pdi + beta.t.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    pdi.trend[t] <- alpha.pdi + beta.t.pdi * t
    mu.pdi[t] <- pdi.trend[t] + ar1.pdi * (PDI[t-1] - pdi.trend[t-1])
  } #t
  
  #Gas Prices
  for(i in 1:2){
    alpha.gas[i] ~ dnorm(0, sd = 10)
    beta.gas[i] ~ dnorm(0, sd = 100)
  }
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  
  #Relative Cost of Gas
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(alpha.gas[era.gas[t]] + beta.gas[era.gas[t]]*t, sd = sig.gas),0,)
    rel.cost[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  } #t
  
  #Unemployment Rate
  sig.une ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
  } #t
  
  #Resident License sales
  sig.res ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(t in 1:n.year){
    res[t] ~ dnorm(0, sd = sig.res)
  }
  
  #Drought Index
  for(r in 1:n.region){
    sig.wpdsi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.pdsi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
      pdsi[t,r] ~ dnorm(0, sd = sig.pdsi[r])
    } #t
  } #r
  
  #Winter Severity
  beta.awssi ~ dnorm(0, sd = 100)
  alpha.awssi ~ dnorm(0, sd = 100)
  for(r in 1:n.region){
    sig.awssi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      for(t in 1:n.year){
      awssi[r,t] ~ dnorm(alpha.awssi + beta.awssi*(era.awssi[t]-1), sd = sig.awssi[r])
    }
  }
  
  # #Rabbits
  # sig.rabbits~ T(dt(0, pow(2.5,-2), 1),0,)
  # beta.rabbits ~ dnorm(0, sd = 2)
  # alpha.rabbits ~ dnorm(0, sd = 1)
  # for(t in 1:n.year){
  #   for(r in 1:n.region){
  #     rabbits[t,r] ~ dnorm(alpha.rabbits + beta.rabbits*t, sd = sig.rabbits)
  #   }
  # }
  
  #Ravens (Highly correlated with Prairie Falcon and RTHawk)
  alpha.rav ~ dnorm(0, 0.001) #intercept
  beta.t.rav ~ dnorm(0, 0.01) #year
  sig.rav~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1.rav ~ dunif(-1,1) #Autoregressive parameter
  
  rav.trend[1] <- alpha.rav + beta.t.rav * 1
  mu.rav[1] <- rav.trend[1]
  for(t in 2:n.year){
    rav.trend[t] <- alpha.rav + beta.t.rav * t
    mu.rav[t] <- rav.trend[t] + ar1.rav * (raven[t-1] - rav.trend[t-1])
  } #t
  
  for(t in 1:n.year){
    raven[t] ~ dnorm(mu.rav[t], sd = sig.rav)
  }
  
  #Northern Harrier
  sig.nhar ~ dunif(0,5)
  for(t in 1:(n.year)){
    nharrier[t] ~ dnorm(0, sd = sig.nhar)
  } #t
} )

data <- list(
  ### Covariates
  une = une, #BL Unemployment information for Nevada, scaled
  res = res, #Residential license sales
  PDI = PDI, #personal disposable income
  GAS = GAS, #gas prices
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  wpdsi = wpdsi, #winter drought index, scaled
  rabbits = rabbits, #Number of rabbits harvested 
  raven = as.vector(bbs.df$raven), #bbs bayes index for ravens, t = 1 is 1975
  nharrier = as.vector(bbs.df$nharrier) #bbs bayes index for northern harriers, t = 1 is 1975
)


### Specify Constants
n.species<- dim(hunters)[1]
sig = rgamma(n.species,1,1)
Lambda = diag(sig)
I = diag(n.species) #identity matrix

sig2 = rgamma(n.species, 1, 1)
Lambda2 = diag(sig2)
I2 = diag(n.species)

constants <- list(
  n.region = 2,
  n.species = nrow(upland),
  n.year = ncol(hunters),
  K = 12,
  
  ### Predictors
  era.gas = c(rep(1,length(1976:2004)),rep(2, length(2005:2017))), #Groupings for change in gas prices 
  era.awssi = c(rep(1,length(1976:1994)),rep(2, length(1995:2001)), rep(1, length(2002:2017))) #Groupings for change in gas prices 
)

### Predictors
une.init <- ifelse(is.na(une), 0, NA)
res.init <- ifelse(is.na(res), 0, NA)
PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), NA)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, NA)
wpdsi.init <- wpdsi
for(i in 1:2){wpdsi.init[,i] <- ifelse(is.na(wpdsi[,i]), 0, NA)}
pdsi.init <- pdsi
for(i in 1:2){pdsi.init[,i] <- ifelse(is.na(pdsi[,i]), 0, NA)}
awssi.init <- awssi
for(i in 1:2){awssi.init[i,] <- ifelse(is.na(awssi[i,]), 0, NA)}
rabbits.init <- rabbits
for(i in 1:2){rabbits.init[,i] <- ifelse(is.na(rabbits[,i]), 0, NA)}
raven.init <- as.vector(ifelse(is.na(bbs.df$raven), 0, NA))
nharrier.init <- as.vector(ifelse(is.na(bbs.df$nharrier), 0, NA))

# Wrapper Function
initsFunction <- function() list(
  ### Predictors
  #Personal Disposable Income
  alpha.pdi = 2.9,
  beta.t.pdi = 0,
  sig.pdi = 1,
  ar1.pdi = 0,
  PDI = PDI.inits,
  #Gas Prices
  GAS = GAS.inits,
  alpha.gas = rep(.9,2),
  beta.gas = rep(0,2),
  sig.gas = 1,
  #Unemployment
  sig.une = 1,
  une = une.init,
  #Resident License Sales
  sig.res = 1,
  res = res.init,
  #Drought Index
  sig.wpdsi = rep(1,2),
  sig.pdsi = rep(1,2),
  wpdsi = wpdsi.init,
  pdsi = pdsi.init,
  #Winter Severity
  sig.awssi = rep(1,2),
  beta.awssi = 0,
  alpha.awssi = 0,
  awssi = awssi.init,
  #Rabbits
  sig.rabbits = 1,
  beta.rabbits = 1,
  alpha.rabbits = 1,
  rabbits = rabbits.init,
  #Ravens
  alpha.rav = 0,
  beta.t.rav = 0,
  sig.rav = 1,
  ar1.rav = 0,
  raven = raven.init,
  #Northern Harrier
  sig.nhar = 1,
  nharrier = nharrier.init
)

inits <- initsFunction()

### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)

model_test$simulate(c())
model_test$initializeInfo()
model_test$calculate()

pars.pred <- c(
    "alpha.pdi",
    "beta.t.pdi",
    "sig.pdi",
    "ar1.pdi",
    "alpha.gas",
    "sig.gas",
    "beta.gas",
    "rel.cost",
    "sig.une",
    "une",
    "sig.res",
    "res",
    "sig.wpdsi",
    "sig.pdsi",
    "sig.awssi",
    "beta.awssi",
    "alpha.awssi",
    "sig.rabbits",
    "beta.rabbits",
    "alpha.rabbits",
    "alpha.rav",
    "beta.t.rav",
    "sig.rav",
    "ar1.rav",
    "sig.nhar"
  )



# Parallel Processing Setup
rm(out)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "pars.pred")) #identify what is to be exported to each cluster

for (j in seq_along(cl)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl[j], "inits")
}

out <- clusterEvalQ(cl, {
  require(nimble)
  require(coda)
  model_test <- nimbleModel( code = code,
                             constants = constants,
                             data =  data,
                             inits = inits )
  
  model_test$simulate(c())
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors =  pars.pred) 
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 200000, niter = 300000, thin = 10)
  
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)

#Find model runtime
end_time <- Sys.time()
end_time - start_time

samples1 <- list(chain1 =  out[[1]], 
                 chain2 =  out[[2]], 
                 chain3 =  out[[3]])

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))

#Save Outputs as file
files <- list(mcmcList1, code)
save(files, file = 'model_output_Predictor.rdata')

# ### Load back
# load(file = 'model_output_Predictor.rdata')
# mcmcList1 <- files[[1]]
# mcmcList2 <- files[[2]]

### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList1, filename = paste("./Traceplots - Predictors Only Model.pdf"))
