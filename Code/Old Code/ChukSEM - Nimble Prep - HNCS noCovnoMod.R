### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2016 #Last year from which data will be used
final.y <- 2017 #Last year to predict
year.hold <- cutoff.y +1

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- final.y - cutoff.y
source("./ChukSEM - Data Prep - Predict.R")


### Run Hunter Effort Solo Model to get estimates of H to create spline inputs
# Sys.time()
# source("./ChukSEM - Hunter Effort Model - Predict.R")
# Sys.time()

load("./model_output_HuntEff_pred.rdata")
mcmcList2 <- files[[1]]

### Load Model Code
source("./ChukSEM - Model Only - HNCS noCovnoMod.R")


### Specify Data Inputs
# Splines
require(splines)
# Function that Constructs B-Spline Base
# from https://github.com/andrewcparnell/jags_examples/blob/master/R%20Code/jags_spline.R
bs_bbase <- function(x, xl = min(x, na.rm = TRUE), xr = max(x, na.rm=TRUE), nseg = 5, deg = 3) {
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

hunter.prime   <- MCMCpstr(mcmcList2, 'H')$H #Extract hunter numbers from Model1

nseg <- 10 #Number of spline segments
BM <- array(NA, dim = c(cut,nseg+3,7,2))
Z  <- array(NA, dim = c(cut,nseg+2,7,2))
D <- diff(diag(ncol(BM[,,1,1])), diff = 1)
Q <- t(D) %*% solve(D %*% t(D))

for(i in 1:7){
  for(j in 1:2){
    BM[,,i,j] <- bs_bbase(hunter.prime[i,,j], nseg = 10)
    Z[,,i,j] <-  BM[,,i,j]%*% Q
  }
}

ZZ <- Z
ZZ[is.na(ZZ)] <- 0

# Time
time <- 1:cut

BM1 <- array(NA, dim = c(cut,nseg+3,7,2))
Z1  <- array(NA, dim = c(cut,nseg+2,7,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))

for(i in 1:7){
  for(j in 1:2){
    BM1[,,i,j] <- bs_bbase(time, nseg = 10)
    Z1[,,i,j] <-  BM1[,,i,j]%*% Q1
  }
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0

data <- list(
  ### Covariates
  # une = une, #BL Unemployment information for Nevada, scaled
  # res = res, #Residential license sales
  # PDI = PDI, #personal disposable income
  # GAS = GAS, #gas prices
  # awssi = awssi, #winter severity index, scaled
  # pdsi = pdsi, #Previous breeding season drought index
  # wpdsi = wpdsi, #winter drought index, scaled
  # rabbits = rabbits, #Number of rabbits harvested 
  # raven = as.vector(bbs.df$raven), #bbs bayes index for ravens, t = 1 is 1975
  # nharrier = as.vector(bbs.df$nharrier), #bbs bayes index for northern harriers, t = 1 is 1975
  
  ### Hunter Effort
  n.hunt = hunters, #Observed number of hunters for each species each year
  Z.hunt = ZZ1, #Spline 
  
  ### Total Harvest
  n.harv = upland,
  Z.harv = ZZ, #Spline
  # 
  # ### Sage Grouse WingBee
  AHY.sg =  wing.b.ahy,
  HY.sg =  wing.b.hy,
  
  ### Chukar Site Abundance
  n.chuk = data.matrix(chukar)
)


### Specify Constants
n.species<- dim(hunters)[1]
sig = rgamma(n.species,1,1)
Lambda = diag(sig)
I = diag(n.species) #identity matrix

sig2 = rgamma(n.species, 1, 1)
Lambda2 = diag(sig2)
I2 = diag(n.species)

constants <- list(
  n.region = 2,
  n.species = nrow(upland),
  n.year = ncol(hunters),
  K = 12,
  
  ### Predictors
  # years.gas = length(1976:2005), #Change point in gas 
  # era.awssi = c(rep(1,length(1976:1994)),rep(2, length(1995:2001)), rep(1, length(2002:2017))), #Groupings for change in gas prices
  
  ### Hunter Effort
  I.hunt = abind(I,I,along = 3),
  
  ### Total Harvest
  I.harv = abind(I2,I2,along = 3),
  
  ### Sage Grouse Wing-Bee
  n.years.sg = n.years.sg,
  rab.use = ifelse(drop.rabbit == "Y", 0, 1),
  
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

## Sage Grouse Wing-Bee
C.sg.i <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.sg.i[,1] <- floor(rowMeans(wing.b.hy, na.rm = T))

hy.sg.init <- wing.b.hy
ahy.sg.init <- wing.b.ahy
for(i in 1:2){
  for(j in 1:ncol(wing.b.ahy)){
    hy.sg.init[i,j] <- ifelse(is.na(wing.b.hy[i,j]), floor(mean(wing.b.hy[i,], na.rm = T)), NA)
    ahy.sg.init[i,j] <- ifelse(is.na(wing.b.ahy[i,j]), floor(mean(wing.b.ahy[i,], na.rm = T)), NA)
  }}
 

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

### Predictors
une.init <- ifelse(is.na(une), 0, NA)
res.init <- ifelse(is.na(res), 0, NA)
PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), NA)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, NA)
wpdsi.init <- wpdsi
for(i in 1:2){wpdsi.init[,i] <- ifelse(is.na(wpdsi[,i]), 0, NA)}
pdsi.init <- pdsi
for(i in 1:2){pdsi.init[,i] <- ifelse(is.na(pdsi[,i]), 0, NA)}
awssi.init <- awssi
for(i in 1:2){awssi.init[i,] <- ifelse(is.na(awssi[i,]), 0, NA)}
rabbits.init <- rabbits
for(i in 1:2){rabbits.init[,i] <- ifelse(is.na(rabbits[,i]), 0, NA)}
raven.init <- as.vector(ifelse(is.na(bbs.df$raven), 0, NA))
nharrier.init <- as.vector(ifelse(is.na(bbs.df$nharrier), 0, NA))

# Wrapper Function
initsFunction <- function() list(
  ### Predictors
  #Personal Disposable Income
  # alpha.pdi = 2.9,
  # beta.t.pdi = 0,
  # sig.pdi = 1,
  # ar1.pdi = 0,
  # # PDI = PDI.inits,
  # #Gas Prices
  # # alpha.gas = rep(.9,2),
  # # beta.t.gas = rep(0,2),
  # sig.gas = 1,
  # mu.gas = .78,
  # # ar1.gas = rep(0,2),
  # GAS = GAS.inits,
  # #Unemployment
  # sig.une = 1,
  # une = une.init,
  # #Resident License Sales
  # sig.res = 1,
  # res = res.init,
  # #Drought Index
  # sig.wpdsi = rep(1,2),
  # sig.pdsi = rep(1,2),
  # wpdsi = wpdsi.init,
  # pdsi = pdsi.init,
  #Winter Severity
  # sig.awssi = rep(1,2),
  # beta.awssi = 0,
  # alpha.awssi = 0,
  # awssi = awssi.init,
  # # #Rabbits
  # # sig.rabbits = 1,
  # # beta.rabbits = 1,
  # # alpha.rabbits = 1,
  # # rabbits = rabbits.init,
  #Ravens
  # alpha.rav = 0,
  # beta.t.rav = 0,
  # sig.rav = 1,
  # ar1.rav = 0,
  # raven = raven.init,
  # #Northern Harrier
  # sig.nhar = 1,
  # nharrier = nharrier.init,

  ###Beta Coefficients
  # mu.drought.hunt = 0,
  # sig.drought.hunt = 1,
  # mu.wintsev.hunt = 0,
  # sig.wintsev.hunt = 1,
  # mu.jobs = 0,
  # sig.jobs = 1,
  # mu.income = 0,
  # sig.income = 1,
  # mu.license = 0,
  # sig.license = 1,
  # 
  # mu.hunters.harv = 0,
  # sig.hunters.harv = 1,
  # mu.wintsev.harv = 0,
  # sig.wintsev.harv = 1,
  # mu.raven = 0,
  # sig.raven = 1,
  # beta.wintsev.harv = rep(0, 7),
  # beta.raven.harv = rep(0, 7),
  # beta.hunters.harv = rep(0, 7),
  beta.spl.harv = array(0, dim = c(7,2,12)),
  # mu.drought.harv = 0,
  # sig.drought.harv = 1,
  # # mu.rabbit = 0,
  # # sig.rabbit = 1,
  # mu.nharrier = 0,
  # sig.nharrier = 1,
  # beta.drought.harv = rep(0, 7),
  # beta.rabbit.harv = rep(0, 7),
  # beta.nharrier.harv = rep(0, 7),
  # ar1.harv = matrix(0, nrow = 7, ncol = 2),
  
  ### Hunter Effort
  n.hunt = n.hunt.i,
  Q.hunt = abind(Q,Q,along = 3),
  # Sigma.hunt = abind(Sigma,Sigma,along = 3),
  P.hunt = abind(P,P,along = 3),
  Lambda.hunt = abind(diag(n.species),diag(n.species),along = 3),
  Delta.hunt = abind(Delta,Delta,along = 3),
  rho.hunt = abind(diag(n.species),diag(n.species),along = 3),
  alpha.hunt = matrix(0, ncol = 2, nrow = 7),
  sig.hunt = matrix(1, ncol = 2, nrow = 7),
  # beta.drought.hunt = rep(0, 7),
  # beta.wintsev.hunt = rep(0, 7),
  # beta.jobs = rep(0, 7),
  # beta.income = rep(0, 7),
  # beta.license = rep(0, 7),
  sig.spl.hunt = matrix(1, ncol = 2, nrow = 7),
  
  ### Total Harvest
  ### Total Harvest
  n.harv = n.harv.i,
  Q.harv = abind(Q2,Q2,along = 3),
  P.harv = abind(P2,P2,along = 3),
  Lambda.harv = abind(diag(n.species),diag(n.species),along = 3),
  Delta.harv = abind(Delta2,Delta2,along = 3),
  rho.harv = abind(diag(n.species),diag(n.species),along = 3),
  alpha.harv = matrix(0, ncol = 2, nrow = 7),
  sig.harv = matrix(1, ncol = 2, nrow = 7),
  log.r.harv = array(0, dim = c(7,(cut)-1,2) ),
  sig.spl.harv = matrix(1, ncol = 2, nrow = 7),
  N = Ni,
  
  ### Sage Grouse Wing-Bee
  theta.sg = rep(1,2),
  # mod.sg = rep(1,2),
  AHY.sg = ahy.sg.init,
  HY.sg = hy.sg.init,
  
  ### Chukar Site Abundance
  theta.chuk = rep(1,2),
  # mod.chuk = rep(1,2),
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
model_test$simulate(c(#'GAS', 'PDI',
                      'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
                      # 'mu.harv.trend',
                      'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',  
                      'theta.sg', 'rate.sg', 'log.r.sg',
                      'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'chuk.eps',#'mod.chuk', 
                      'BPH'))
model_test$initializeInfo()
model_test$calculate()


pars1 <- c(### Hunter Effort
  "alpha.hunt",
  # "beta.drought.hunt",
  # "beta.wintsev.hunt",
  # "beta.jobs",
  # "beta.income",
  # "beta.license",
  "pred.spl.hunt",
  
  ### Total Harvest
  "alpha.harv",
  # "beta.drought.harv",
  # "beta.wintsev.harv",
  # "beta.hunters.harv",
  # "beta.raven.harv",
  # "beta.nharrier.harv",
  # "ar1.harv",
  "pred.spl.harv",
  
  ### Sage Grouse Wing-Bee
  # "alpha.sg",
  # "beta.wintsev.sg",
  # "beta.drought.sg",
  # "beta.rabbit.sg",
  # "beta.raven.sg",
  # "beta.nharrier.sg",
  # "mod.sg",
  "theta.sg",
  
  ### Chukar Site Abundance
  # "alpha.chuk",
  # "beta.wintsev.chuk",
  # "beta.drought.chuk",
  # "beta.rabbit.chuk",
  # "beta.raven.chuk",
  # # "beta.nharrier.chuk",
  # "mod.chuk",
  "theta.chuk"
)

pars2 <- c(### Hunter Effort
  "H",
  # "rho.hunt",
  "log.r.hunt",
  
  ### Total Harvest
  "N",
  # "rho.harv",
  "log.r.harv",
  
  # ### Sage Grouse Wing-Bee
  # "log.r.sg",
  # 
  # ### Chukar Site Abundance
  # 
  # "log.r.chuk",
  
  ### Birds per Hunter
  "BPH")

pars.pred <- c(
  # "alpha.pdi",
  # "beta.t.pdi",
  # "sig.pdi",
  # "ar1.pdi",
  # "alpha.gas",
  # "beta.t.gas",
  # "mu.gas",
  # "sig.gas",
  # "ar1.gas",
  # "rel.cost",
  # "sig.une",
  # "une",
  # "sig.res",
  # "res",
  # "sig.wpdsi",
  # "sig.pdsi",
  # "sig.awssi",
  # "beta.awssi",
  # "alpha.awssi",
  # "sig.rabbits",
  # "beta.rabbits",
  # "alpha.rabbits",
  # "alpha.rav",
  # "beta.t.rav",
  # "sig.rav",
  # "ar1.rav",
  # "sig.nhar",
  # "raven",
  # "awssi"
  
)

# Parallel Processing Setup
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
  
  model_test$simulate(c(#'GAS', 'PDI',
    'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
    # 'mu.harv.trend',
    'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',  
    'theta.sg', 'rate.sg', 'log.r.sg',
    'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
    'BPH'))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors =  c(pars1, pars2))#, monitors2 = pars.pred)
  # mcmcConf <-  configureMCMC( model_test,   monitors =  pars1, monitors2 = pars2)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 95000, niter = 100000, thin = 5, thin2 = 5)
  # samplesList <- runMCMC(Cmcmc,nburnin = 5, niter = 100, thin = 1, thin2 = 1)
  
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)

#Find model runtime
end_time <- Sys.time()
end_time - start_time

samples1 <- list(chain1 =  out.full.predict[[1]],
                 chain2 =  out.full.predict[[2]],
                 chain3 =  out.full.predict[[3]])

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))

# samples1 <- list(chain1 =  out.full.predict[[1]]$samples,
#                  chain2 =  out.full.predict[[2]]$samples,
#                  chain3 =  out.full.predict[[3]]$samples)
# 
# mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
# 
# samples2 <- list(chain1 =  out.full.predict[[1]]$samples2,
#                  chain2 =  out.full.predict[[2]]$samples2,
#                  chain3 =  out.full.predict[[3]]$samples2)
# 
# mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

#Save Outputs as file
files <- list(mcmcList1, code)
# files <- list(mcmcList1, mcmcList2, code)
save(files, file = paste("./Holdout ", year.hold, '/model_output_FullModel_predict.rdata', sep = ""))

# 
# mcmcList1 <- files[[1]]
# mcmcList2 <- files[[2]]



source("./ChukSEM - Estimate Check.R")


### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList1, filename = paste("./Holdout ", year.hold, '/TraceOut - Full.pdf', sep = ""))

# MCMCtrace(mcmcList2, filename = paste("./Holdout ", year.hold, '/TraceOut - Full - Predictors.pdf', sep = ""))
# 
# 


