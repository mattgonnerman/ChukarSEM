### Specify Data Inputs
require(splines)
bs_bbase <- function(x, xl = min(x, na.rm = TRUE), xr = max(x, na.rm=TRUE), nseg = 10, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(bs(x, knots = knots, degree = deg, Boundary.knots = c(knots[1], knots[length(knots)])), nrow = length(x))
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  return(bs_matrix)
}

nseg <- 15 #Number of spline segments
time <- 1:cut
B <- bs_bbase(x = time, nseg = nseg)

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
  bobcat = bobcat.df[,2],
 
  ### Harvest Data
  n.hunt = hunters/1000, #Observed number of hunters for each species each year
  basis = B, #Spline Base
  n.harv = upland/1000,
  # bph.survey = surveybph,
  
  ### Chukar Site Abundance
  n.chuk = data.matrix(chukar)
)


### Specify Constants
n.species<- dim(hunters)[1]
sig = rgamma(n.species,1,1)
Lambda <- Lambda2 <- diag(sig)
I2 <- I <- diag(n.species) #identity matrix

constants <- list(
  n.region = 2,
  n.species = nrow(upland),
  n.year = ncol(hunters),
  K = dim(B)[2],

  ### Predictors
  era.awssi = c(rep(0,length(1975:1994)),rep(1, length(1995:2001)), rep(0, length(2002:final.y))), #Groupings for change in gas prices

  ### Hunter Effort
  I.hunt = abind(I,I,along = 3),

  ### Total Harvest
  I.harv = abind(I2,I2,along = 3),
  
  ### BPH Survey
  # n.cut = length(1976:cutoff.y),
  
  ### Chukar Site Abundance
  n.site = nrow(chukar),
  n.year.chuk = ncol(chukar),
  reg.chuk = chuk.reg
)


### Specify Initial Values
## Hunter Effort
nu = n.species + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.species, n.species)
for (j in 1:n.species){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)

n.hunt.i <- ifelse(is.na(hunters), floor(mean(hunters, na.rm = T)), NA)


##Total Harvest
nu = n.species + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.species, n.species)
for (j in 1:n.species){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)

nu = n.species + 1
Q2 = rinvwishart(nu, I)    # note output is an array
Delta2 = matrix(0, n.species, n.species)

for (j in 1:n.species){
  Delta2[j,j] = Q2[j,j]^(-0.5)
}
P2 = Delta2 %*% Q2 %*% Delta2
Sigma2 = Lambda2 %*% P2 %*% Lambda2

n.hunt.i <- ifelse(is.na(hunters), floor(mean(hunters, na.rm = T)), NA)
n.harv.i <- ifelse(is.na(upland), floor(mean(upland, na.rm = T)), NA)

Ni <- array(NA, c(nrow(n.harv.i), ncol(n.harv.i), 2))
Ni[,1,] <- upland[,1,] + 50

## Chukar Site Abundance
chukar_na <- chukar

for(i in 1:nrow(chukar_na)){
  for(j in 1:ncol(chukar_na)){
    if(is.na(chukar_na[i,j])){
      chukar_na[i,j] <- floor(mean(as.matrix(chukar[i,]), na.rm = T))
    }else{
      chukar_na[i,j] <- NA
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
C.chuk.init[,2:ncol(C.chuk.init)] <- NA

econ.inits <- econ_data %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))
bbs.inits <- as.data.frame(bbs.df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))
awssi.inits <- as.matrix(as.data.frame(awssi.df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))
pdsi.inits <- as.matrix(as.data.frame(pdsi_df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))
bobcat.inits <- as.matrix(as.data.frame(bobcat.df[,2]) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))

# Wrapper Function
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
  bobcat = as.numeric(bobcat.inits),
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
  sig.bob = 1,
  
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
  mu.bobcat = 0,
  sig.bobcat = 1,
  beta.bobcat.harv = rep(0, n.species),
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
  
  # sig.bph = rep(.1, n.species),
  
  ### Chukar Site Abundance
  theta.chuk = rep(1,2),
  mod.chuk = rep(1,2),
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
model_test$simulate(c(
  'mu.hunt', 'beta.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
  'mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',
  'theta.chuk','rate.chuk', ' log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
  'BPH', 'BPH2'))
model_test$initializeInfo()
model_test$calculate()


estimates <- c(
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
  "theta.chuk",
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
  'beta.awssi',
  
  'sig.drought',
  
  'sig.bob'
)

# Parallel Processing Setup
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
    'mu.hunt', 'beta.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
    'mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',
    'theta.chuk','rate.chuk', ' log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
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
MCMCtrace(mcmcList1, filename = "./Output/TraceOut - Estimates.pdf")
MCMCtrace(mcmcList2, filename = "./Output/TraceOut - Covariates.pdf")

