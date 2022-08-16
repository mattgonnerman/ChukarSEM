#####################################################################################
### Data
data <- list(
  ### Covariates
  awssi = awssi.df, #winter severity index, scaled
  ravens = bbs.df[,1],
  rthawk = bbs.df[,2],
  nharr = bbs.df[,3],
  pfal = bbs.df[,4],
  une = econ_data[,2], #BL Unemployment information for Nevada, scaled
  res = econ_data[,3], #Resident Licenses
  pdi = econ_data[,4], #Personal Disposable Income
  gas = econ_data[,1], #Gas Prices in May
  pdsi = pdsi_df,
  
  ### Hunter Effort
  n.hunt = chuk_hunt/1000, #Observed number of hunters for each species each year
  basis = B, #Spline Base
  I.hunt = abind(I,I,along = 3),
  
  ### Total Harvest
  n.harv = chuk_harv/1000,
  I.harv = abind(I2,I2,along = 3),
  
  ### Chukar Site Abundance
  n.chuk = data.matrix(survey.abun)
)

#####################################################################################
### Constants
constants <- list(
  n.region = 2,
  n.counties = n.counties,
  n.year = ncol(chuk_hunt),
  K = dim(B)[2],
  
  ### Chukar Site Abundance
  n.site = nrow(chukar),
  n.year.chuk = ncol(chukar),
  reg.chuk = chuk.reg
)

#####################################################################################
### Initial Values
initsFunction <- function() list(
  ### Covariates
  gas = econ.inits$Gas.May,
  une = econ.inits$Une,
  res = econ.inits$Licenses,
  pdi = econ.inits$PDI,
  ravens = bbs.inits$raven,
  rthawk = bbs.inits$rthawk,
  nharr = bbs.inits$nharrier,
  pfal = bbs.inits$pfalcon,
  awssi = awssi.inits,
  pdsi = pdsi.inits,
  sig.econ.pred = rep(1,4),
  mu.econ.p = rep(0,4),
  zeta.econ = rep(1,ncol(econ_data)),
  pred.econ.prime = rep(0, cut),
  mu.econ.pred = matrix(0, 4, cut),
  sig.bbs.pred = rep(1,4),
  mu.bbs.p = rep(0,4),
  zeta.bbs = rep(1,ncol(bbs.df)),
  pred.bbs.prime = rep(0, cut),
  mu.bbs.pred = matrix(0, 4, cut),
  beta.awssi = 0,
  alpha.awssi = 0,
  sig.awssi = rep(1,2),
  sig.drought = rep(1,2),
  
  ### Covariates
  mu.econ = 0,
  sig.econ = 1,
  beta.econ.hunt = rep(0, n.species),
  beta.spl.hunt = array(0, dim = c(n.species,2,dim(B)[2])),
  sig.spl.hunt = matrix(1, ncol = 2, nrow = n.species),
  
  mu.wintsev.harv = 0,
  sig.wintsev.harv = 1,
  beta.wintsev.harv = rep(0, n.species),
  mu.bbs = 0,
  sig.bbs = 1,
  beta.bbs.harv = rep(0, n.species),
  mu.pdsi = 0,
  sig.pdsi = 1,
  beta.pdsi.harv = rep(0, n.species),
  mu.hunter.harv = c(0,0),
  sig.hunter.harv = c(1,1),
  beta.hunter.harv = matrix(0, n.species, 2),
  
  ### Hunter Effort
  sig.H =  matrix(1, n.species,2),
  n.hunt = n.hunt.i/1000,
  Q.hunt = abind(Q,Q,along = 3),
  P.hunt = abind(P,P,along = 3),
  Lambda.hunt = abind(diag(n.species),diag(n.species),along = 3),
  Delta.hunt = abind(Delta,Delta,along = 3),
  rho.hunt = abind(diag(n.species),diag(n.species),along = 3),
  alpha.hunt = matrix(0, ncol = 2, nrow = n.species),
  sig.hunt = matrix(1, ncol = 2, nrow = n.species),
  
  ### Total Harvest
  sig.N =  matrix(1, n.species ,2),
  n.harv = n.harv.i/1000,
  Q.harv = abind(Q2,Q2,along = 3),
  P.harv = abind(P2,P2,along = 3),
  Lambda.harv = abind(diag(n.species),diag(n.species),along = 3),
  Delta.harv = abind(Delta2,Delta2,along = 3),
  rho.harv = abind(diag(n.species),diag(n.species),along = 3),
  alpha.harv = matrix(0, ncol = 2, nrow = n.species),
  sig.harv = matrix(1, ncol = 2, nrow = n.species),
  log.r.harv = array(0, dim = c(n.species,(cut)-1,2) ),
  N = Ni,
  
  ### Chukar Site Abundance
  theta.chuk = rep(1,2),
  mod.chuk = rep(1,2),
  mu.chuk = c(0,0),
  n.chuk = as.matrix(chukar_na),
  log.r.chuk = r.chuk.init
)

inits <- initsFunction()


#####################################################################################
### Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)
model_test$simulate(c(
  'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
  'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',
  'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
  'BPH', 'BPH2'
))
model_test$initializeInfo()
model_test$calculate()


#####################################################################################
### Parameters to Montior
pars1 <- c(### Hunter Effort
  "alpha.hunt",
  "beta.econ.hunt",
  "beta.spl.hunt",
  "Sigma.hunt",
  
  ### Total Harvest
  "alpha.harv",
  "beta.wintsev.harv",
  "beta.bbs.harv",
  "beta.pdsi.harv",
  # "beta.spl.harv",
  "Sigma.harv"
)

pars2 <- c(
  ### Hunter Effort
  "H",'sig.H',
  "rho.hunt",
  "log.r.hunt",
  ### Total Harvest
  "N",'sig.N',
  "log.r.harv",
  'rho.harv',
  # ### Chukar Site Abundance
  "log.r.chuk",
  'mod.chuk',
  ### Birds per Hunter
  "BPH",'BPH2'
)

pars.pred <- c(
  'pred.econ.prime',
  'pred.bbs.prime',
  'sig.econ.pred', 
  'sig.bbs.pred',
  
  'mu.bbs.p',
  'zeta.econ',
  'zeta.bbs',
  'alpha.awssi', 
  'beta.awssi',
  'awssi',
  
  'beta.econ.hunt',
  'sig.H',
  'alpha.hunt',
  'beta.spl.hunt',
  'latent.trend',
  
  'beta.wintsev.harv',
  'beta.bbs.harv',
  'beta.pdsi.harv',
  'beta.hunter.harv',
  'sig.N','mu.harv',
  'mu.hunt',
  'mod.chuk'
  
)


#####################################################################################
### Run Model (Parallel Processing)
rm(out)
rm(out.full.predict)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "pars1", "pars2", "pars.pred")) #identify what is to be exported to each cluster

for (j in seq_along(cl)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl[j], "inits")
}

out.full.predict <- clusterEvalQ(cl, {
  require(nimble)
  require(coda)
  
  model_test <- nimbleModel( code = code,
                             constants = constants,
                             data =  data,
                             inits = inits )
  model_test$simulate(c(
    'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
    'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',
    'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
    'BPH', 'BPH2'
  ))
  model_test$initializeInfo()
  model_test$calculate()
  
  mcmcConf <-  configureMCMC( model_test,   monitors =  pars2, monitors2 = c(pars1, pars.pred))
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  samplesList <- runMCMC(Cmcmc,nburnin = 50000, niter = 125000, thin = 5, thin2 = 5)
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time


#####################################################################################
### Save Outputs
samples1 <- list(chain1 =  out.full.predict[[1]]$samples,
                 chain2 =  out.full.predict[[2]]$samples,
                 chain3 =  out.full.predict[[3]]$samples)

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))

samples2 <- list(chain1 =  out.full.predict[[1]]$samples2,
                 chain2 =  out.full.predict[[2]]$samples2,
                 chain3 =  out.full.predict[[3]]$samples2)

mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

#Save Outputs as file
files <- list(mcmcList1, mcmcList2, code, data)
save(files, file = "./Output/NDOW_Upland_SEM_output.rdata")

#Traceplots
MCMCtrace(mcmcList1, filename = "./Output/TraceOut - Full.pdf")
MCMCtrace(mcmcList2, filename = "./Output/TraceOut - Full - Predictors.pdf")

