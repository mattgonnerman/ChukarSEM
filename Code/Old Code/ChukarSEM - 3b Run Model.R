# Parameters monitored
pars1 <- c('b0')

start_time <- Sys.time()

inits <- initsFunction()

#Check Model before running fully


library(parallel)
library(coda)
set.seed(2)
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
  
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  c('lbo1',
                                                           'lbo',
                                                           'N', 
                                                           'lambda',
                                                           'pred',
                                                           'beta.drought2',
                                                           'beta.drought',
                                                           'beta.jobs',
                                                           'hunt.eps', 
                                                           'H',
                                                           'bph', #birds per hunter
                                                           'mu',
                                                           'mu2',
                                                           'pred1',
                                                           'rho',
                                                           'rho2',
                                                           'lambda1', 
                                                           'beta.pressure',
                                                           'beta.trend',
                                                           'lbo2', 
                                                           'chuk.eps', 
                                                           'lambda2',
                                                           'C', 
                                                           'lambda2',
                                                           'theta2',
                                                           'log.r2', 
                                                           'sigma.chuk'))
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  # samplesList <- runMCMC(Cmcmc,nburnin = 250000, niter = 500000, thin = 10, thin2 = 10)
  samplesList <- runMCMC(Cmcmc,nburnin = 2500, niter = 5000, thin = 10, thin2 = 10)
  
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

library(coda)
library(MCMCvis)

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))



mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'ChukarSEM_model_output.rdata')
