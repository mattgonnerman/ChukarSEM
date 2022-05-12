### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Model
# Based on Boyce et al. 2001? https://www.jstor.org/stable/3803103
# https://georgederpa.github.io/teaching/countModels.html
code <- nimbleCode( {
  ################################################################################
  ### Sage Grouse Wing-Bee ###
  beta.drought.sg ~ dnorm(0, 0.01)
  beta.wintsev.sg ~ dnorm(0, 0.01)
  beta.rabbit.sg ~ dnorm(0, 0.01)
  beta.raven.sg ~ dnorm(0, 0.01)
  beta.nharrier.sg ~ dnorm(0, 0.01)
  
  for(r in 1:n.region){
    alpha.sg[r] ~ dnorm(0, sd = 100) #Intercept
    # mod.sg[r] ~ dlogis(0,1) #Constant modifier to translate harv change to wingb change
    theta.sg[r] ~ dunif(0,1) #NB "probability" parameter, between 0 and 1
    
    for(t in 1:n.years.sg){
      # sg.eps[r, t] <- mod.sg[r] * log.r.harv[7, t+28, r] #log.r.harv[t=25] is 2004
      
      log.r.sg[r,t] <- alpha.sg[r] + 
        # sg.eps[r,t] + #Unlinked change in recruitment
        beta.drought.sg * pdsi[t+28,r] + #previous breeding season drought index
        beta.wintsev.sg * awssi[r,t+28] + #previous winter severity
        beta.rabbit.sg * rabbits[t+27,r] + #previous year's rabbit harvest
        beta.raven.sg * raven[t+29] + #previous spring BBS raven index
        beta.nharrier.sg * nharrier[t+29]  #previous spring BBS northern harrier index
      
      rate.sg[r,t] <- theta.sg[r]/(theta.sg[r] + (AHY.sg[r,t]*exp(log.r.sg[r,t]))) #NB rate
      HY.sg[r,t] ~ dnegbin(size = rate.sg[r,t], prob = theta.sg[r])
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
  n.years.sg = n.years.sg,
  n.region = 2
)


### Specify Initial Values
C.sg.i <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.sg.i[,1] <- floor(rowMeans(wing.b.hy, na.rm = T))

initsFunction <- function() list(   
  ### Sage Grouse Wing-Bee
  theta.sg = rep(1,2),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  mod.sg = rep(1,2),
  alpha.sg =  rep(0,2),
  beta.drought.sg = 0,
  beta.wintsev.sg = 0,
  beta.rabbit.sg = 0,
  beta.raven.sg = 0,
  beta.nharrier.sg = 0
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

#Parallel Processing Setup
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
start_time <- Sys.time()
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000)

clusterExport(cl, c("code", "inits", "data", "constants", "pars1"))

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
  
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  pars1)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
 samplesList <- runMCMC(Cmcmc,nburnin = 40000, niter = 60000, thin = 3, thin2 = 3)
  
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
MCMCtrace(mcmcList2, filename = "Traceplots - SG WingBee MCMC.pdf")