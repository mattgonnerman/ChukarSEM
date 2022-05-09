### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Model
# Based on Boyce et al. 2001? https://www.jstor.org/stable/3803103
# https://georgederpa.github.io/teaching/countModels.html
code <- nimbleCode( {
  ################################################################################
  ### Sage Grouse Wing-Bee ###
  for(r in 1:n.region){
    mod.sg[r] ~ dlogis(0,1) #Constant modifier to translate harv change to wingb change
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    sigma.chuk[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #Amount of between year variation in change
    
    for(t in 1:n.years.sg){
      log.r.sg[r,t] ~ dnorm(log.r.sg[7, t+24, r]*mod.sg[r], sd = sigma.chuk[r])
      lambda.sg[r,t] <- exp(log.r.sg[r,t]) #Change in wing-b is proportional to change in harvest
      C.HY.sg[r,t] <- dpois(AHY.sg[r,t]*lambda.sg[r,t])
      rate.sg[r,t] <- theta.sg[r]/(theta.sg[r] + C.HY.sg[r,t]) #NB rate
      HY.sg[r,t] ~ dnegbin(rate.sg[r,t], theta.sg[r])
      lambda.sg[r,t-1] <- exp(log.r.sg[r,t-1]) # Tranform between finite and instantaneous growth rate
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
  rabbits = rabbits #Number of rabbits harvested 
             )


### Contansts
constants <- list(
  n.years.sg = n.years.sg,
  n.region = 2
)


### Specify Initial Values
C.sg.i <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.sg.i[,1] <- floor(rowMeans(wing.b, na.rm = T))

initsFunction <- function() list(   
  ### Sage Grouse Wing-Bee
  theta.sg = rep(1,2),
  # sigma.sg = rep(0,2),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  C.sg = C.sg.i,
  # mod.sg = rep(1,2),
  alpha.sg =  rep(0,2),
  beta.drought.sg = 0,
  beta.wintsev.sg = 0,
  beta.rabbit.sg = 0
  # alpha.sg =  rep(0,2),
  # beta.drought.sg = rep(0, 2),
  # beta.wintsev.sg = rep(0, 2),
  # beta.rabbit.sg = rep(0, 2)
)

inits <- initsFunction()



### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           inits = inits,
                           constants = constants,
                           data =  data)

model_test$simulate(c('log.r.sg',
                      'C.sg'))
model_test$initializeInfo()
model_test$calculate()

#Set Monitors
pars1 <- c("alpha.sg",
           "beta.wintsev.sg",
           "beta.drought.sg",
           "beta.rabbit.sg",
           
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
  
 samplesList <- runMCMC(Cmcmc,nburnin = 40000, niter = 60000, thin = 10, thin2 = 10)
  
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