### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Model
code <- nimbleCode( {
  ################################################################################
  ### Chukar Site Abundance ###
  beta.drought.chuk ~ dnorm(0, 0.01)
  beta.wintsev.chuk ~ dnorm(0, 0.01)
  beta.rabbit.chuk ~ dnorm(0, 0.01)
  
  for(p in 1:n.site){
    alpha.chuk[p] ~ dnorm(0, sd = 100) #Intercept
    
    theta.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    # sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
    # mod.chuk[p] ~ dlogis(0,1)
    
    for(t in 2:n.year.chuk){
      # mu.chuk[p, t-1] <- mod.chuk[p] * log.r.harv[3, t+12, reg.chuk[p]] #log.r.harv[t=14] is 1990-1991
      # chuk.eps[p,t-1]  ~ dnorm(mu.chuk[p, t-1], sd = sigma.chuk[p]) #Coeffient for annual change
      
      log.r.chuk[p,t-1] <- alpha.chuk[p] + 
        # chuk.eps[p,t-1] + #Unlinked change in abundance, log.harv.chuk[t=1] is 1990-1991
        beta.drought.chuk * pdsi[t+13,chuk.reg[p]] + #previous breeding season drought index
        beta.wintsev.chuk * awssi[chuk.reg[p],t+14] + #previous year's winter severity
        beta.rabbit.chuk * rabbits[t+13,chuk.reg[p]] #previous year's rabbit harvest
      lambda.chuk[p,t-1] <- exp(log.r.chuk[p,t-1]) #Change in abundance
      
      C.chuk[p,t] <- lambda.chuk[p,t-1] * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[p]/(theta.chuk[p] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(rate.chuk[p,t-1], theta.chuk[p]) #obs. # of chukars follow neg-bin
    } #t
  } #p 
})


### Specify Data Inputs
data <- list(
  ### Chukar Site Abundance
  n.chuk = data.matrix(chukar),
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  rabbits = rabbits #Number of rabbits harvested 
)



### Specify Contansts
constants <- list(
  n.site = 13,
  n.year.chuk = ncol(chukar),
  chuk.reg = chuk.reg
)


### Initial Values
chukar_na <- chukar 
chukar_na <- ifelse(is.na(chukar == TRUE), floor(mean(as.matrix(chukar), na.rm = T)), NA)

initsFunction <- function() list( 
  ### Chukar Site Abundance
  theta.chuk = rep(1,13),
  # sigma.chuk = rep(0,2),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.chuk-1),
  n.chuk = chukar_na,
  # mod.chuk = rep(1,2),
  alpha.chuk =  rep(0,13),
  beta.drought.chuk = 0,
  beta.wintsev.chuk = 0,
  beta.rabbit.chuk = 0
  # alpha.sg =  rep(0,2),
  # beta.drought.sg = rep(0, 2),
  # beta.wintsev.sg = rep(0, 2),
  # beta.rabbit.sg = rep(0, 2)
)

inits <- initsFunction()


### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)

model_test$simulate(c("log.r.chuk", "C.chuk", "rate.chuk"))
model_test$initializeInfo()
model_test$calculate()

### Monitors
pars1 <- c(
  "alpha.chuk",
  "beta.wintsev.chuk",
  "beta.drought.chuk",
  "beta.rabbit.chuk",
  
  "theta.chuk",
  "log.r.chuk"
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



mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

#Save Outputs as file
files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'model_output_ChukSite.rdata')

### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList2, filename = "Traceplots - Chuk Site MCMC.pdf")