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
  # pdsi = pdsi_df,
  # bobcat = bobcat.df[,2],
  
  ### Hunter Effort
  n.hunt = chuk_hunt/1000, #Observed number of hunters for each species each year
  basis = B, #Spline Base
  I.hunt = I,
  
  ### Total Harvest
  n.harv = chuk_harv/1000,
  I.harv = I2,
  
  # bph.survey = surveybph,
  
  ### Chukar Site Abundance
  n.chuk = data.matrix(survey.abun)
)

#####################################################################################
### Constants
constants <- list(
  #Harvest Information
  n.region = 2,
  n.counties = n_county,
  n.year = ncol(chuk_hunt),
  K = dim(B)[2],
  reg.county = county_reg,
  # n.cut = length(1990:cutoff.y),
  
  ### Chukar Site Abundance
  n.site = nrow(chuk_hunt),
  n.year.chuk = ncol(survey.abun),
  county.site = survey.county,
  
  ### Predictors
  era.awssi = c(rep(0,length(1990:1994)),rep(1, length(1995:2001)), rep(0, length(2002:final.y))) #Groupings for change in gas prices
  
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
  # pdsi = pdsi.inits,
  # bobcat = bobcat.inits,
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
  # sig.drought = rep(1,2),
  # sig.bob = 1,
  
  ### Covariates
  mu.econ = 0,
  sig.econ = 1,
  beta.econ.hunt = rep(0, n_county),
  beta.spl.hunt = matrix(0, nrow = n_county, ncol = dim(B)[2]),
  sig.spl.hunt = rep(1, n_county),
  
  sig.wintsev.harv = 1,
  beta.wintsev.harv = 0,
  sig.bbs = 1,
  beta.bbs.harv = 0,
  mu.hunter.harv = 0,
  sig.hunter.harv = 1,
  beta.hunter.harv = rep(0, n_county),

  ### Hunter Effort
  sig.H =  rep(1, n_county),
  n.hunt = n.hunt.i/1000,
  Q.hunt = Q,
  P.hunt = P,
  Lambda.hunt = diag(n_county),
  Delta.hunt = Delta,
  rho.hunt = diag(n_county),
  alpha.hunt = rep(0, n_county),
  sig.hunt = rep(1, n_county),
  
  ### Total Harvest
  sig.N =  rep(1, n_county),
  n.harv = n.harv.i/1000,
  Q.harv = Q2,
  P.harv = P2,
  Lambda.harv = diag(n_county),
  Delta.harv = Delta2,
  rho.harv = diag(n_county),
  alpha.harv = rep(0, n_county),
  sig.harv = rep(1, n_county),
  log.r.harv = matrix(0, nrow = n_county, ncol = (cut)-1),
  N = Ni,
  
  # sig.bph = .1,
  
  ### Chukar Site Abundance
  theta.chuk = 1,
  mod.chuk = rep(1,n_county),
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
estimates <- c(
  ### Hunter Effort
  "H",
  "rho.hunt",
  ### Total Harvest
  "N",'sig.N',
  "log.r.harv",
  'rho.harv',
  # ### Chukar Site Abundance
  "log.r.chuk",
  "theta.chuk",
  'n.chuk',
  ### Birds per Hunter
  "BPH",'BPH2'
)

coefficients <- c(
  ### Hunter Effort
  "alpha.hunt",
  "beta.econ.hunt",
  "beta.spl.hunt",
  "latent.trend",
  
  ### Total Harvest
  "alpha.harv",
  "beta.wintsev.harv",
  "beta.hunter.harv",
  # "beta.bobcat.harv",
  "beta.bbs.harv",
  # "beta.pdsi.harv",
  
  ### Chukar Site Abundance
  'theta.chuk',
  # 'mu.site.chuk',
  'mod.chuk'
)

predictors <- c(
  'pred.econ.prime',
  'pred.bbs.prime',
  'mu.econ.p',
  'mu.bbs.p',
  'sig.econ.pred', 
  'sig.bbs.pred',
  'zeta.econ',
  'zeta.bbs',
  
  'alpha.awssi', 
  'beta.awssi'
  
  # 'sig.drought',
  
  # 'sig.bob'
)


#####################################################################################
### Run Model (Parallel Processing)
rm(out.full.predict)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "estimates", "coefficients", "predictors")) #identify what is to be exported to each cluster

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
  
  mcmcConf <-  configureMCMC( model_test,   monitors =  estimates, monitors2 = c(coefficients, predictors))
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  samplesList <- runMCMC(Cmcmc,nburnin = 50000, niter = 125000, thin = 5, thin2 = 5)
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)

#Find model runtime
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
files <- list(mcmcList1, mcmcList2, code, data, county_order, county_reg)
save(files, file = "./Output/NDOW_ChukOnly_SEM_output.rdata")

#Traceplots
MCMCtrace(mcmcList1, filename = "./Output/ChukOnly - Estimates.pdf")

MCMCtrace(mcmcList2, filename = "./Output/ChukOnly - Covariates.pdf")

