### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2017 #Only need to change this to adjust the number of years

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- 2017-cutoff.y
source("./ChukSEM - Data Prep - Predict.R")

### Model
# Based on Boyce et al. 2001? https://www.jstor.org/stable/3803103
# https://georgederpa.github.io/teaching/countModels.html
code <- nimbleCode( {
  ################################################################################
  ### Predictors
  #Drought Index
  for(r in 1:n.region){
    # sig.wpdsi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.pdsi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(t in 1:n.year){
      # wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
      pdsi[t,r] ~ dnorm(0, sd = sig.pdsi[r])
    } #t
  } #r
  
  #Winter Severity
  sig.awssi~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.awssi[1] ~ dnorm(0, 0.01)
  beta.awssi[2] ~ dnorm(0, 0.01)
  for(r in 1:n.region){
    alpha.awssi[r] ~ dnorm(1.5, 1)
    for(t in 1:n.year){
      awssi[r,t] ~ dnorm(alpha.awssi[r] + beta.awssi[era.awssi[t]]*t, sd = sig.awssi)
    }
  }
  
  #Rabbits
  sig.rabbits~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.rabbits ~ dnorm(0, 0.01)
  for(r in 1:n.region){
    alpha.rabbits[r] ~ dnorm(0, 1)
    for(t in 1:n.year){
      rabbits[t,r] ~ dnorm(alpha.rabbits[r] + beta.rabbits*t, sd = sig.rabbits)
    }
  }
  
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
  
  ### Sage Grouse Wing-Bee ###
  beta.drought.sg ~ dnorm(0, 0.01)
  beta.wintsev.sg ~ dnorm(0, 0.01)
  beta.rabbit.sg ~ dnorm(0, 0.01)
  beta.raven.sg ~ dnorm(0, 0.01)
  beta.nharrier.sg ~ dnorm(0, 0.01)
  
  for(r in 1:n.region){
    alpha.sg[r] ~ dnorm(0, sd = 100) #Intercept
    # mod.sg[r] ~ dlogis(0,1) #Constant modifier to translate harv change to wingb change
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    
    for(t in 1:n.years.sg){
      # sg.eps[r, t] <- mod.sg[r] * log.r.harv[7, t+28, r] #log.r.harv[t=25] is 2004
      
      log.r.sg[r,t] <- alpha.sg[r] + 
        # sg.eps[r,t] + #Unlinked change in recruitment
        beta.drought.sg * pdsi[t+28,r] + #previous breeding season drought index
        beta.wintsev.sg * awssi[r,t+27] + #previous winter severity (2003-2004 winter for t = 1)
        beta.rabbit.sg * rabbits[t+27,r] + #previous year's rabbit harvest
        beta.raven.sg * raven[t+29] + #previous spring BBS raven index
        beta.nharrier.sg * nharrier[t+29]  #previous spring BBS northern harrier index
      
      rate.sg[r,t] <- theta.sg[r]/(theta.sg[r] + (AHY.sg[r,t]*exp(log.r.sg[r,t]))) #NB rate
      HY.sg[r,t] ~ dnegbin(prob = rate.sg[r,t], size = theta.sg[r])
    } #t
  } #r
})


### Specify Data Inputs
data <- list(
  ### Sage Grouse WingBee
  AHY.sg =  wing.b.ahy,
  HY.sg =  wing.b.hy,
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  rabbits = rabbits, #Number of rabbits harvested 
  raven = as.vector(bbs.df$raven), #bbs bayes index for ravens, t = 1 is 1975
  nharrier = as.vector(bbs.df$nharrier) #bbs bayes index for northern harriers, t = 1 is 1975
             )


### Contansts
constants <- list(
  n.year = ncol(hunters),
  n.years.sg = n.years.sg,
  n.region = 2,
  era.awssi = c(rep(1,length(1976:1994)),rep(2, length(1995:2001)), rep(1, length(2002:2017))) #Groupings for change in gas prices 
  
)


### Specify Initial Values
C.sg.i <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.sg.i[,1] <- floor(rowMeans(wing.b.hy, na.rm = T))

hy.sg.init <- wing.b.hy
ahy.sg.init <- wing.b.ahy
for(i in 1:2){
  for(j in 1:ncol(wing.b.ahy)){
    hy.sg.init[i,j] <- ifelse(is.na(wing.b.hy[i,j]), floor(mean(wing.b.hy[i,], na.rm = T)), NA)
    ahy.sg.init[i,j] <- ifelse(is.na(wing.b.ahy[i,j]), floor(mean(wing.b.ahy[i,], na.rm = T)), NA)
  }}

### Predictors
une.init <- ifelse(is.na(une), 0, NA)
res.init <- ifelse(is.na(res), 0, NA)

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

initsFunction <- function() list(   
  #Drought Index
  sig.pdsi = rep(1,2),
  pdsi = pdsi.init,
  #Winter Severity
  sig.awssi = 1,
  beta.awssi = rep(0,2),
  alpha.awssi = rep(0,2),
  awssi = awssi.init,
  #Rabbits
  sig.rabbits = 1,
  beta.rabbits = 0,
  alpha.rabbits = rep(0,2),
  rabbits = rabbits.init,
  #Ravens
  alpha.rav = 0,
  beta.t.rav = 0,
  sig.rav = 1,
  ar1.rav = 0,
  raven = raven.init,
  #Northern Harrier
  sig.nhar = 1,
  nharrier = nharrier.init,
  ### Sage Grouse Wing-Bee
  theta.sg = rep(1,2),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  # mod.sg = rep(1,2),
  alpha.sg =  rep(0,2),
  beta.drought.sg = 0,
  beta.wintsev.sg = 0,
  beta.rabbit.sg = 0,
  beta.raven.sg = 0,
  beta.nharrier.sg = 0,
  AHY.sg = ahy.sg.init,
  HY.sg = hy.sg.init
)

inits <- initsFunction()



### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           inits = inits,
                           constants = constants,
                           data =  data)

model_test$simulate(c('theta.sg', 'rate.sg', 'log.r.sg'))
model_test$initializeInfo()
model_test$calculate()

#Set Monitors
pars1 <- c("alpha.sg",
           "beta.wintsev.sg",
           "beta.drought.sg",
           "beta.rabbit.sg",
           "beta.raven.sg",
           "beta.nharrier.sg",
           
           "theta.sg",
           "log.r.sg"
)

pars2 <- c(#"alpha.pdi",
           # "beta.t.pdi",
           # "sig.pdi",
           # "ar1.pdi",
           # "alpha.gas",
           # "sig.gas",
           # "beta.gas",
           # "rel.cost",
           # "sig.une",
           # "une",
           # "sig.res",
           # "res",
           # "sig.wpdsi",
           "sig.pdsi",
           # "wdpsi",
           "pdsi", 
           "sig.awssi",
           "beta.awssi",
           "alpha.awssi",
           "awssi",
           "sig.rabbits",
           "beta.rabbits",
           "alpha.rabbits",
           "rabbits",
           "alpha.rav",
           "beta.t.rav",
           "sig.rav",
           "ar1.rav",
           "raven",
           "sig.nhar",
           "nharrier")

#Parallel Processing Setup
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
start_time <- Sys.time()
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000)

clusterExport(cl, c("code", "inits", "data", "constants", "pars1", "pars2"))

for (j in seq_along(cl)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl[j], "inits")
}

out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  model_test <- nimbleModel( code = code, 
                             constants = constants,  
                             data =  data, 
                             inits = inits )
  
  mcmcConf <-  configureMCMC( model_test,   monitors =  pars1, monitors2 = pars2)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
 samplesList <- runMCMC(Cmcmc,nburnin = 40000, niter = 60000, thin = 5, thin2 = 5)
  
  return(samplesList)
})

end_time <- Sys.time()
end_time - start_time

stopCluster(cl)

samples2 <- list(chain1 =  out[[1]]$samples2, 
                 chain2 =  out[[2]]$samples2, 
                 chain3 =  out[[3]]$samples2)

samples1    <- list(chain1 =  out[[1]]$samples, 
                    chain2 =  out[[2]]$samples, 
                    chain3 =  out[[3]]$samples)



mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

#Save Outputs as file
files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'model_output_SGWingBee.rdata')

### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots

MCMCtrace(mcmcList1, filename = "Traceplots - SG WingBee MCMC - Process.pdf")
MCMCtrace(mcmcList2, filename = "Traceplots - SG WingBee MCMC - Predictors.pdf")
