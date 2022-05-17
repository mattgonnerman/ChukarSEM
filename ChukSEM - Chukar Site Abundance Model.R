### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Model
code <- nimbleCode( {
  ################################################################################
  ### Chukar Site Abundance ###
  beta.drought.chuk ~ dnorm(0, 0.01)
  beta.wintsev.chuk ~ dnorm(0, 0.01)
  beta.rabbit.chuk ~ dnorm(0, 0.01)
  beta.raven.chuk ~ dnorm(0, 0.01)
  beta.nharrier.chuk ~ dnorm(0, 0.01)
  
  for(r in 1:n.region){
    alpha.chuk[r] ~ dnorm(0, sd = 100) #Intercept
    theta.chuk[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    # mod.chuk[r] ~ dlogis(0,1)
  }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.chuk){
      # chuk.eps[p, t-1] <- mod.chuk[chuk.reg[p]] * log.r.harv[3, t+13, reg.chuk[p]] #log.r.harv[t=13] is 1990-1991
      
      log.r.chuk[p,t-1] <- alpha.chuk[chuk.reg[p]] + 
        # chuk.eps[p,t-1] + #Unlinked change in abundance, log.harv.chuk[t=1] is 1990-1991
        beta.drought.chuk * pdsi[t+12,chuk.reg[p]] + #previous breeding season drought index
        beta.wintsev.chuk * awssi[chuk.reg[p],t+12] + #previous year's winter severity
        beta.rabbit.chuk * rabbits[t+12,chuk.reg[p]] + #previous year's rabbit harvest
        beta.raven.chuk * raven[t+13] + #previous year's spring BBS index
        beta.nharrier.chuk * nharrier[t+13] #previous year's spring BBS index
      
      C.chuk[p,t] <- exp(log.r.chuk[p,t-1]) * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[chuk.reg[p]]/(theta.chuk[chuk.reg[p]] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(prob = rate.chuk[p,t-1], size = theta.chuk[chuk.reg[p]]) #obs. # of chukars follow neg-bin
    } #t
  } #p 
})


### Specify Data Inputs
data <- list(
  ### Chukar Site Abundance
  n.chuk = data.matrix(chukar),
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  rabbits = rabbits, #Number of rabbits harvested 
  raven = as.vector(bbs.df$raven), #bbs bayes index for ravens, t = 1 is 1975
  nharrier = as.vector(bbs.df$nharrier) #bbs bayes index for northern harriers, t = 1 is 1975
)



### Specify Contansts
constants <- list(
  n.site = 13,
  n.region = 2,
  n.year.chuk = ncol(chukar),
  chuk.reg = chuk.reg
)


### Initial Values
chukar_na <- chukar 

for(i in 1:nrow(chukar_na)){
  for(j in 1:ncol(chukar_na)){
    if(is.na(chukar_na[i,j])){
      chukar_na[i,j] <- floor(mean(as.matrix(chukar[i,]), na.rm = T))
    # }else{
    #   chukar_na[i,j] <- NA
    }
  }
}

r.chuk.init <- matrix(NA, nrow = nrow(chukar), ncol = ncol(chukar)-1)
for(i in 1:nrow(chukar_na)){
  for(j in 2:ncol(chukar_na)){
    r.chuk.init[i,j-1] <- chukar_na[i,j]/chukar_na[i,j-1]
  }
}

C.chuk.init <- chukar
C.chuk.init[is.na(C.chuk.init)] <- floor(mean(as.matrix(chukar), na.rm = T))

initsFunction <- function() list( 
  ### Chukar Site Abundance
  # theta.chuk = rep(1,13),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.chuk-1),
  # mod.chuk = rep(1,2),
  alpha.chuk =  rep(0,2),
  beta.drought.chuk = 0,
  beta.wintsev.chuk = 0,
  beta.rabbit.chuk = 0,
  beta.raven.chuk = 0,
  beta.nharrier.chuk = 0,
  n.chuk = as.matrix(chukar_na),
  log.r.chuk = r.chuk.init
)

inits <- initsFunction()


### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)

model_test$simulate(c('C.chuk', 'rate.chuk', 'theta.chuk'))
model_test$initializeInfo()
model_test$calculate()

### Monitors
pars1 <- c(
  "alpha.chuk",
  "beta.wintsev.chuk",
  "beta.drought.chuk",
  "beta.rabbit.chuk",
  "beta.raven.chuk",
  "beta.nharrier.chuk",
  
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