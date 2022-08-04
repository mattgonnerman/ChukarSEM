###################################
### Sage Grouse Wing-Bee Model ###
##################################
set.seed(2)
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "coda", "MCMCvis", "parallel"),
       require, character.only = T)

### Model
code <- nimbleCode( {
  ################################################################################
  ### Sage Grouse Wing-Bee ###
  ### Chukar Site Abundance
  for(r in 1:n.region){
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.sg[r,1] ~ dpois(wing.b[r,1]) #Equivalent of Poisson lambda
    
    alpha.sg[r] ~ dlogis(0,1) #Intercept for
    sigma.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod.sg[r] ~ dlogis(0,1)
    
    for(t in 2:n.years.sg){
      mu.sg[r, t-1] <- mod.sg[r]
      sg.eps[r,t-1] ~ dnorm(mu.sg[r, t-1], sd = sigma.sg[r])
      
      log.r.sg[r,t-1] <- alpha.sg[r] + sg.eps[r,t-1] #Variation in growth
      lambda.sg[r,t-1] <- exp(log.r.sg[r,t-1]) # Tranform between finite and instantaneous growth rate
      
      C.sg[r,t] <- lambda.sg[r,t-1] * C.sg[r,t-1] #Change in Count over time
      
      rate.sg[r,t-1] <- theta.sg[r]/(theta.sg[r] + C.sg[r,t]) #NB rate
      wing.b[r,t] ~ dnegbin(rate.sg[r,t-1], theta.sg[r]) #Wing-Bee data follows negative binomial
    } #t
  } #r
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
C.sg.i <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.sg.i[,1] <- floor(rowMeans(wing.b, na.rm = T))

initsFunction <- function() list(   
  ### Sage Grouse Wing-Bee
  alpha.sg =  rep(0,2),
  theta.sg = rep(1,2),
  sigma.sg = rep(0,2),
  sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  C.sg = C.sg.i,
  mod.sg = rep(1,2)
)

inits <- initsFunction()

### Contansts
constants <- list(
  n.years.sg = n.years.sg,
  n.region = 2
)

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

#Parallel Processing Setup
start_time <- Sys.time()

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
  
 samplesList <- runMCMC(Cmcmc,nburnin = 5000, niter = 10000, thin = 10, thin2 = 10)
  
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
