lapply(c("parallel", "coda"), require, character.only = T)

set.seed(2)

start_time <- Sys.time() # To track runtime

nc <- 3    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "pars1")) #identify what is to be exported to each cluster

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
  
  model_test$simulate(c('sig', 'pred1', 'mu', 'une', 'wpdsi',  'beta.trend'))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  c('lbo1', #regression intercept?
                                                           'beta.trend', #Spline smoothing factor
                                                           'beta.general', #general license sales
                                                           'beta.drought2', #Change in hunters associated with drought
                                                           'beta.jobs', #Change in hunters associated with employment
                                                           'beta.income',
                                                           'wpdsi', #drought conditions
                                                           'une', #unemployment
                                                           'GAS', #gas price
                                                           'PDI', #personal disposable income
                                                           'REL.COST', #GAS/PDI
                                                           'hunt.eps', # lambda for poisson describing Hunter numbers(unlinked regression output)
                                                           'H', #log(hunt.eps)
                                                           'mu', #mean lambda for poisson dscribing hunter numbers 
                                                           'rho', #correlation coefficent
                                                           'pred1', #Spline smoothing 
                                                           'ar1', #change in PDI over time?
                                                           'lambda1')) #Change in hunter numbers
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 100, niter = 1000, thin = 1, thin2 = 1)
  
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)

#Find model runtime
end_time <- Sys.time()
end_time - start_time

samples2 <- list(chain1 =  out[[1]]$samples2, 
                 chain2 =  out[[2]]$samples2, 
                 chain3 =  out[[3]]$samples2)

samples1    <- list(chain1 =  out[[1]]$samples, 
                    chain2 =  out[[2]]$samples, 
                    chain3 =  out[[3]]$samples)

