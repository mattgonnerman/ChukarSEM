### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Total Harvest ###
  for(r in 1:n.region){
    mu.drought.harv[r] ~ dnorm(0, 0.01)
    sig.drought.harv[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.wintsev.harv[r] ~ dnorm(0, 0.01)
    sig.wintsev.harv[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.rabbit[r] ~ dnorm(0, 0.01)
    sig.rabbit[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    
    for(s in 1:n.species){
      # Regression Coefficients
      alpha.harv[s,r] ~ dnorm(0, sd = 1)
      beta.drought.harv[s,r] ~ dnorm(mu.drought.harv[r], sd = sig.drought.harv[r])
      beta.wintsev.harv[s,r] ~ dnorm(mu.wintsev.harv[r], sd  = sig.wintsev.harv[r])
      beta.rabbit[s,r] ~ dnorm(mu.rabbit[r], sd  = sig.rabbit[r])
      for(k in 1:K){
        beta.spl.harv[s,r,k] ~ dnorm(0, sd = sig.spl.harv[s,r])
      } #k
      sig.spl.harv[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Total harvest, Year 1
      
      for(t in 2:(n.year)){
        mu.harv[s,t-1,r] <- alpha.harv[s,r] +#regression formula
          beta.drought.harv[s,r] * pdsi[t,r] + #previous breeding season drought index
          beta.wintsev.harv[s,r] * awssi[r,t] + #concurrent winter severity
          beta.rabbit[s,r] * rabbits[t,r] + #concurrent rabbit harvest
          inprod(beta.spl.harv[s,r,1:K], Z.harv[t,1:K,s,r]) #spline smoothing

        pred.spl.harv[s,r,t-1] <- inprod(beta.spl.harv[s,r,1:K], Z.harv[t,1:K,s,r]) #Derive spline smoothing for examination later
        
        lambda.harv[s,t-1,r] <- exp(log.r.harv[s,t-1,r]) #link function
        N[s,t,r] <- lambda.harv[s,t-1,r] * N[s,t-1,r] #number available = change since last year
        n.harv[s,t,r] ~  dpois(N[s,t,r]) #Number harvested follows Poisson
      } #t
    } #s 
    
    for(t in 2:(n.year)){
      #Change in total harvest, log.r.harv[t=1] is 1976-1977
      log.r.harv[1:n.species,t-1,r]  ~ dmnorm(mu.harv[1:n.species,t-1,r],
                                              cov =  Sigma.harv[1:n.species,1:n.species,r])
    } #t
    
    ### Correlation Matrices
    Q.harv[1:n.species,1:n.species,r] ~ dinvwish(S = I.harv[1:n.species,1:n.species,r], df = n.species + 1)
    
    for(s in 1:n.species){
      sig.harv[s,r] ~ dgamma(1,1)
      Delta.harv[s,s,r] <- pow(Q.harv[s,s,r], -0.5)
      Lambda.harv[s,s,r] <- sig.harv[s,r]
    } #s
    
    for (s1 in 2:n.species){
      for (s2 in 1:(s1-1)){
        Lambda.harv[s1,s2,r] <- 0
        Delta.harv[s1,s2,r] <- 0
      } #s2
    } #s1
    
    Sigma.harv[1:n.species,1:n.species,r] <- Lambda.harv[1:n.species,1:n.species,r] %*% P.harv[1:n.species,1:n.species,r] %*% Lambda.harv[1:n.species,1:n.species,r]  
    P.harv[1:n.species,1:n.species,r] <- Delta.harv[1:n.species,1:n.species,r] %*% Q.harv[1:n.species,1:n.species,r] %*% Delta.harv[1:n.species,1:n.species,r]  
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho.harv[s1,s2,r] <- Sigma.harv[s1,s2,r]/sqrt(Sigma.harv[s1,s1,r] * Sigma.harv[s2,s2,r])   
      } #s2
    } #s1
  } #r
  
})

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

nseg <- 10
BM1 <- array(NA, dim = c(cut,nseg+3,7,2))
Z1  <- array(NA, dim = c(cut,nseg+2,7,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))
time <- 1:(cut)
for(i in 1:7){
  for(j in 1:2){
    BM1[,,i,j] <- bs_bbase(time, nseg = 10)
    Z1[,,i,j] <-  BM1[,,i,j]%*% Q1
    
  }
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0

data <- list(
  ### Total Harvest
  n.harv = upland,
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  rabbits = rabbits, #Number of rabbits harvested 
  Z.harv = ZZ1 #Spline
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
  n.species = 7,
  n.year = ncol(hunters),
  K = 12,
  
  ### Total Harvest
  I.harv = abind(I2,I2,along = 3)
  
)


### Specify Initial Values
## Hunter Effort/Total Harvest
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
n.chuk.i <- ifelse(is.na(chukar), floor(mean(as.matrix(chukar), na.rm = T)), NA)
chukar_na <- chukar 
chukar_na <- ifelse(is.na(chukar == TRUE), 1, 0)

Xi <- chukar 
for(i in 1:ncol(chukar_na)){
  for(j in 1:nrow(chukar_na)){
    if(chukar_na[j,i] == 1){
      chukar_na[j,i] <- floor(rnorm(1, rowMeans( chukar[j,], na.rm = TRUE), 5))
    }
  }
}
for(j in 1:nrow(chukar_na)){
  if(is.na(Xi[j,1]) == TRUE){
    chukar_na[j,1] <- chukar_na[j,1]
  } else
    chukar_na[j,1] <- Xi[j,1]
}

chukar_na[chukar_na ==0] <- NA

## Sage Grouse Wing-Bee
C.sg.i <- matrix(NA, nrow = 2, ncol = n.years.sg)
C.sg.i[,1] <- floor(rowMeans(wing.b, na.rm = T))

# Wrapper Function
initsFunction <- function() list(   
  ### Predictors
  sig.spl.harv = matrix(1, ncol = 2, nrow = 7),
  mu.drought.harv = rep(0,2),
  sig.drought.harv = rep(1,2),
  mu.wintsev.harv = rep(0,2),
  sig.wintsev.harv = rep(1,2),
  mu.rabbit = rep(0,2),
  sig.rabbit = rep(1,2),
  
  ### Total Harvest
  n.harv = n.harv.i,
  Q.harv = abind(Q2,Q2,along = 3),
  P.harv = abind(P2,P2,along = 3),
  Lambda.harv = abind(diag(n.species),diag(n.species),along = 3),
  Delta.harv = abind(Delta2,Delta2,along = 3),
  rho.harv = abind(diag(n.species),diag(n.species),along = 3),
  alpha.harv = matrix(0, ncol = 2, nrow = 7),
  sig.harv = matrix(1, ncol = 2, nrow = 7),
  log.r.harv = array(0, dim = c(7,cut-1,2) ),
  N = Ni,
  beta.drought.harv = matrix(0, 7, 2),
  beta.wintsev.harv = matrix(0, 7, 2),
  beta.rabbit = matrix(0, 7, 2),
  beta.spl.harv = array(0, dim = c(7,2,12))
)

inits <- initsFunction()


### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c("mu.harv"))
model_test$initializeInfo()
model_test$calculate()

#Set Monitors
# pars1 <- c(
#   ### Total Harvest
#   "mu.drought.harv",
#   "sig.drought.harv",
#   "mu.wintsev.harv",
#   "sig.wintsev.harv",
#   "mu.rabbit",
#   "sig.rabbit",
#   
#   "alpha.harv",
#   "beta.drought.harv",
#   "beta.wintsev.harv",
#   "beta.spl.harv",
#   "beta.rabbit",
#   "sig.spl.harv",
#   
#   "pred.spl.harv",
#   
#   "log.r.harv",
#   
#   "Q.harv",
#   "sig.harv",
#   "rho.harv"
# )
pars1 <- c(
  "alpha.harv",
  "beta.drought.harv",
  "beta.wintsev.harv",
  "beta.spl.harv",
  "beta.rabbit"
)

### Parallel Processing Code
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
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
  
  model_test$simulate(c("alpha.sg", "sg.eps", "C.sg", "theta.sg", "mod.sg", "rate.sg"))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  pars1)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 40000, niter = 60000, thin = 10, thin2 = 10)
  
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

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

#Save Outputs as file
files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'model_output_TotHarv.rdata')

### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList2, filename = "Traceplots - Total Harvest MCMC.pdf")