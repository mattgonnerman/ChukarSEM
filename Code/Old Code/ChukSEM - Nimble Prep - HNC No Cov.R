### Load Packages and Set Data Subset Info
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2016 #Last year from which data will be used
cutoff.y.chuk <- cutoff.y + 1 #How many years of Chukar site data do we have (default to just the spring of same year)
final.y <- 2017 #Last year to predict
year.hold <- cutoff.y +1
drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models
n.add.y <- final.y - cutoff.y

### Run Hunter Effort Solo Model to get estimates of H to create spline inputs
source("E:/GitHub/ChukarSEM/ChukSEM - Data Prep.R")
source("./ChukSEM - Hunter Effort Model - Predict.R")


###########################
### Only Run Above Once ###
###########################


### Load Model Code
source("./ChukSEM - Model Only - HNC No Cov.R")


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

hunter.prime   <- MCMCpstr(mcmcList1, 'H')$H #Extract hunter numbers from Model1

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
  ### Hunter Effort
  n.hunt = hunters, #Observed number of hunters for each species each year
  Z.hunt = ZZ1, #Spline 
  
  ### Total Harvest
  n.harv = upland,
  Z.harv = ZZ, #Spline
  
  ### Chukar Site Abundance
  n.chuk = data.matrix(chukar)
)


### Specify Constants
constants <- list(
  n.region = 2,
  n.species = nrow(upland),
  n.year = ncol(hunters),
  K = 12,
  
  ### Hunter Effort
  I.hunt = abind(I,I,along = 3),
  
  ### Total Harvest
  I.harv = abind(I2,I2,along = 3),
  
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

# Wrapper Function
initsFunction <- function() list(
  ### Hunter Effort
  n.hunt = n.hunt.i,
  Q.hunt = abind(Q,Q,along = 3),
  P.hunt = abind(P,P,along = 3),
  Lambda.hunt = abind(diag(n.species),diag(n.species),along = 3),
  Delta.hunt = abind(Delta,Delta,along = 3),
  rho.hunt = abind(diag(n.species),diag(n.species),along = 3),
  alpha.hunt = matrix(0, ncol = 2, nrow = 7),
  sig.hunt = matrix(1, ncol = 2, nrow = 7),
  sig.spl.hunt = matrix(1, ncol = 2, nrow = 7),
  
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
                      'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt',
                      'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',  
                      'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
                      'BPH'
                      ))
model_test$initializeInfo()
model_test$calculate()


pars1 <- c(### Hunter Effort
  "alpha.hunt",
  
  ### Total Harvest
  "alpha.harv",
  
  ### Chukar Site Abundance
  "mod.chuk",
  "theta.chuk"
)

pars2 <- c(### Hunter Effort
  "H",
  "rho.hunt",
  "log.r.hunt",
  
  ### Total Harvest
  "N",
  "rho.harv",
  "log.r.harv",
  
  ### Birds per Hunter
  "BPH")

# Parallel Processing Setup
rm(out)
rm(out.full.predict)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "pars1", "pars2")) #identify what is to be exported to each cluster

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
    'BPH'
  ))
  model_test$initializeInfo()
  model_test$calculate()
  
  mcmcConf <-  configureMCMC( model_test,   monitors =  c(pars1, pars2))
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  samplesList <- runMCMC(Cmcmc,nburnin = 95000, niter = 100000, thin = 5, thin2 = 5)
  
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

#Save Outputs as file
files <- list(mcmcList1, code)
save(files, file = paste("./Holdout ", year.hold, '/model_output_FullModel_predict.rdata', sep = ""))

#Traceplots
MCMCtrace(mcmcList1, filename = paste("./Holdout ", year.hold, '/TraceOut - Full.pdf', sep = ""))

#Preliminary Graphs
source("./ChukSEM - Estimate Check.R")