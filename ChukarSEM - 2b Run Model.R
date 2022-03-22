lapply(c("parallel", "coda"), require, character.only = T)

### Check Model Code
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c('sig',
                      'pred1',
                      'mu',
                      'REL.COST',
                      'mu.pdi',
                      'une',
                      'wpdsi',
                      'beta.trend'))
model_test$initializeInfo()
model_test$calculate()


### Parallel Processing Code
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
                                                           'beta.income', #Change in relative cost of gas
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
  
  samplesList <- runMCMC(Cmcmc,nburnin = 4000, niter = 6000, thin = 1, thin2 = 1)
  
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

### Check that model is running correctly
# This creates the uncompiled model based on the model likelihood (code), the data constants
# (constants), the primary data (data) and the initial values (inits)
model_test <- nimbleModel(code = code, constants = constants,  data =  data, inits = inits )

# Once the initial model is built with no errors, you can use the simulation function to 
# simulate new initial values for all parameters (estimates or derived) within the model. 
# Using simulate with c('parameter1', 'parameter2') will simulate just the specified 
# parameters, where simulate() will simulate new values for all parameters. The latter can 
# cause you issues as it overwrites the specified initial values from the inits list 
# provided, which may be required for certain models to run. The primary error you may run
# into here is if there is a dimension mismatch between the likelihood and initial values, 
# which needs to be corrected for the model to run.
model_test$simulate(c('sig', 'pred1', 'mu', 'une', 'wpdsi', 'beta.trend'))

# This code alerts us to which parameters (derived and estimated) are currently missing 
# initial values. Although it is possible for a model to run with no issues when initial 
# values are not provided for the estimated parameters (and the model run is even less 
# sensitive to derived parameters), we will need all parameters to have initial values to 
# complete the last step in this 'model validation' step. I simply keep adding the parameters
# listed in this step to the simulation list above until I get an error (and fix it) or 
# confirmation that all parameters have initial values (this step can be a headache).
model_test$initializeInfo()

# Lastly, this step attempts to calculate the -log likelihood of the entire model based on 
# its current initial values (no different than running a model for a single iteration with 
# no effort to 'update' the model). If an -Inf is returned, the -log likelihood can't be 
# calculated because there the model/data combination is not valid, this error indicates that
# some part of the model will never estimate (generally, the parameter effected (and any 
# parameter that 'touches' the problematic parameter will never accept a new value during 
# the Bayesian updating process, so the initial values will be returned.

# If an NA is returned, it indicates that it couldn't calculate the -log likelihood
# because the initial values were incomplete, this provides us with no information 
# whether the model/data combination is valid or not, hence the important of simulate/
# initalizeInfo step. This model could be valid or not.

# If a negative number is returned, that's the -logLikelihood value from the model, 
# which indicates that the model was able to successfully calculate it based on the 
# initial conditions provided. This is what we want to see.
model_test$calculate()
