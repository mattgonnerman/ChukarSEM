###################################
### Sage Grouse Wing-Bee Model ###
##################################
set.seed(2)
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "coda", "MCMCvis", "parallel"),
       require, character.only = T)

### Model
code <- nimbleCode( {
  ### Chukar Site Abundance
  for(r in 1:n.reg){
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.sg[r,1] ~ dpois(wing.b[r,1]) #Equivalent of Poisson lambda
    
    lbo.sg[r] ~ dlogis(0,1) #Intercept for 
    sigma.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod.sg[r] ~ dlogis(0,1)
    
    for(t in 1:(n.years.sg-1)){
      # sg.eps[r,t-1] <- dnorm(mod.sg[r] * log.r[3, t + ts.sg, 1], sd = sigma.sg[r]) 
      sg.eps[r,t] ~ dnorm(mod.sg[r], sd = sigma.sg[r]) 
      log.r.sg[r,t] <- lbo.sg[r] + sg.eps[r,t] #Variation in growth
      lambda.sg[r,t] <- exp(log.r.sg[r,t]) # Tranform between finite and instantaneous growth rate
    
      C.sg[r,t+1] <- lambda.sg[r,t] * C.sg[r,t] #Change in Count over time
      
      rate.sg[r,t] <- theta.sg[r]/(theta.sg[r] + C.sg[r,t+1])
      wing.b[r,t+1] ~ dnegbin(rate.sg[r,t], theta.sg[r])
      
    } 
  }
})



### Data Preparation ###
#Format WingBee Data
sg.wingb <- read.csv("./Data/SG_WingData_2004-2020.csv") %>%
  select(Region = NDOWREGION, Year, AHY.Male, AHY.Female, HY.Male, HY.Female) %>%
  mutate(Total = rowSums(.[3:6], na.rm = T)) %>%
  select(Region, Year, Total) %>%
  group_by(Region, Year) %>%
  summarize(Total.SG = sum(Total)) %>%
  pivot_wider(names_from = "Region", values_from = "Total.SG")
  
wing.b <- t(sg.wingb[,-1])
n.years.sg <- ncol(wing.b)
time.shift.sg <- min(sg.wingb$Year) - min(harvest_data$Year) - 1

# Wrap in List
data <- list(wing.b =  wing.b
)

### Initial Values
C.init <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.init[] <- wing.b[] + 50

rate.sg.init <- 1/(1+wing.b[,1:(n.years.sg-1)])

sg.eps.init <- matrix(0, nrow = 2, ncol = n.years.sg-1)

inits <- list( 
  C.sg = C.init,
  log.r.sg =  matrix(1, nrow = 2, ncol = n.years.sg - 1),
  theta.sg = rep(1,2),
  lbo.sg =  rep(0,1),
  mod.sg = rep(1,2),
  sigma.sg = rep(0,2),
  rate.sg = rate.sg.init,
  sg.eps = sg.eps.init
)

### Contansts
constants <- list(
  n.years.sg = n.years.sg,
  ts.sg = time.shift.sg,
  n.reg = 2
)

### Monitors
pars1 <- c(
  'lbo.sg', 
  'sg.eps',
  'C.sg', 
  'lambda.sg',
  'theta.sg',
  'log.r.sg', 
  'sigma.sg'
)

### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           inits = inits,
                           constants = constants,
                           data =  data)

model_test$simulate(pars1)
model_test$initializeInfo()
model_test$calculate()

#Parallel Processing Setup
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
  
 samplesList <- runMCMC(Cmcmc,nburnin = 500, niter = 10000, thin = 10, thin2 = 10)
  
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

files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'ChukarSEM_model_output.rdata')