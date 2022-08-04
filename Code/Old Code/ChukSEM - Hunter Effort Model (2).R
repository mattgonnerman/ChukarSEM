### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Economic Predictors ###
  for(i in 1:4){
    sig.econ.pred[i] ~ dgamma(1,1)
    mu.econ.p[i] ~ dnorm(0, 1)
  }
  zeta.econ[1] <- 1
  for(i in 2:4){
    zeta.econ[i] ~ dnorm(0, 1)
  }
  for(t in 1:n.year){
    pred.econ.prime[t] ~ dnorm(0, 1) # Latent Predator Index
    gas[t] ~ dnorm(mu.econ.pred[1,t], sd = sig.econ.pred[1])
    une[t] ~ dnorm(mu.econ.pred[2,t], sd = sig.econ.pred[2])
    res[t] ~ dnorm(mu.econ.pred[3,t], sd = sig.econ.pred[3])
    pdi[t] ~ dnorm(mu.econ.pred[4,t], sd = sig.econ.pred[4])
    
    mu.econ.pred[1,t] <- mu.econ.p[1] + zeta.econ[1] * pred.econ.prime[t]
    mu.econ.pred[2,t] <- mu.econ.p[2] + zeta.econ[2] * pred.econ.prime[t]
    mu.econ.pred[3,t] <- mu.econ.p[3] + zeta.econ[3] * pred.econ.prime[t]
    mu.econ.pred[4,t] <- mu.econ.p[4] + zeta.econ[4] * pred.econ.prime[t]
  }
  
  ################################################################################
  ### Hunter Effort ###
  mu.econ ~ dnorm(0, sd = 100)
  sig.econ ~ T(dt(0, pow(2.5, -2), 1), 0, )
  for (s in 1:n.species) {
    beta.econ.hunt[s] ~ dnorm(mu.econ, sd = sig.econ)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.hunt[s,r] ~ dnorm(5, sd = 12) #sd = 1)
      for(k in 1:K){
        beta.spl.hunt[s,r,k] ~ dnorm(0, sd = sig.spl.hunt[s,r])
      } #k
      sig.spl.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(t in 1:n.year){ 
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,t,r] <- alpha.hunt[s,r] + #intercept
          beta.econ.hunt[s] * pred.econ.prime[t] + #SEM economic indicator
          inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #spline smoothing
        
        H[s,t,r] <- exp(hunt.eps[s,t,r]) #Log Link
        
        n.hunt[s,t,r] ~ dpois(H[s,t,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.hunt[s,t,r] <- H[s,t+1,r]/H[s,t,r]
        log.r.hunt[s,t,r] <- log(lambda.hunt[s,t,r])
      }
    } #s
    
    for(t in 1:n.year){
      hunt.eps[1:n.species,t,r] ~ dmnorm(mu.hunt[1:n.species,t,r], cov =  Sigma.hunt[1:n.species,1:n.species,r] )
    }
    
    # Correlation Matrices
    Q.hunt[1:n.species,1:n.species,r] ~ dinvwish(S = I.hunt[1:n.species,1:n.species,r], df = n.species + 1)
    
    for(s in 1:n.species){
      sig.hunt[s,r] ~ dgamma(1,1)
      Delta.hunt[s,s,r] <- pow(Q.hunt[s,s,r], -0.5)
      Lambda.hunt[s,s,r] <- sig.hunt[s,r]
    } #s
    
    for (s1 in 2:n.species){
      for (s2 in 1:(s1-1)){
        Lambda.hunt[s1,s2,r] <- 0
        Delta.hunt[s1,s2,r] <- 0
      } #s2
    } #s1
    
    Sigma.hunt[1:n.species,1:n.species,r] <- Lambda.hunt[1:n.species,1:n.species,r] %*% P.hunt[1:n.species,1:n.species,r] %*% Lambda.hunt[1:n.species,1:n.species,r]  
    P.hunt[1:n.species,1:n.species,r] <- Delta.hunt[1:n.species,1:n.species,r] %*% Q.hunt[1:n.species,1:n.species,r] %*% Delta.hunt[1:n.species,1:n.species,r]
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho.hunt[s1,s2,r] <- Sigma.hunt[s1,s2,r]/sqrt(Sigma.hunt[s1,s1,r] * Sigma.hunt[s2,s2,r])   
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

nseg <- 10 #Number of spline segments
# Time
time <- 1:cut

BM1 <- array(NA, dim = c(cut,nseg+3,5,2))
Z1  <- array(NA, dim = c(cut,nseg+2,5,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))

for(i in 1:5){
  for(j in 1:2){
    BM1[,,i,j] <- bs_bbase(time, nseg = 10)
    Z1[,,i,j] <-  BM1[,,i,j]%*% Q1
  }
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0

data <- list(
  ### Covariates
  une = econ_data[,2], #BL Unemployment information for Nevada, scaled
  res = econ_data[,3], #Resident Licenses
  pdi = econ_data[,4], #Personal Disposable Income
  gas = econ_data[,1], #Gas Prices in May
  
  ### Hunter Effort
  n.hunt = hunters, #Observed number of hunters for each species each year
  Z.hunt = ZZ1 #Spline
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

  ### Hunter Effort
  I.hunt = abind(I,I,along = 3)
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

### Predictors
econ.inits <- econ_data %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))

# Wrapper Function
initsFunction <- function() list(
  ### Predictors
  sig.econ.pred = rep(1,4),
  mu.econ.p = rep(0,4),
  zeta.econ = rep(1,ncol(econ_data)),
  gas = econ.inits$Gas.May,
  une = econ.inits$Une,
  res = econ.inits$Licenses,
  pdi = econ.inits$PDI,
  pred.econ.prime = rep(0, cut),
  mu.econ.pred = matrix(0, 4, cut),
  
  ###Beta Coefficient Means/SDs
  mu.econ = 0,
  sig.econ = 1,
  beta.econ.hunt = rep(0, 5),
  
  ### Hunter Effort
  n.hunt = n.hunt.i,
  Q.hunt = abind(Q,Q,along = 3),
  P.hunt = abind(P,P,along = 3),
  Lambda.hunt = abind(diag(n.species),diag(n.species),along = 3),
  Delta.hunt = abind(Delta,Delta,along = 3),
  rho.hunt = abind(diag(n.species),diag(n.species),along = 3),
  alpha.hunt = matrix(0, ncol = 2, nrow = 5),
  sig.hunt = matrix(1, ncol = 2, nrow = 5),
  sig.spl.hunt = matrix(1, ncol = 2, nrow = 5)
)

inits <- initsFunction()

### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)
model_test$simulate(c(
  'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'mu.econ.pred',
  'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt'
))
model_test$initializeInfo()
model_test$calculate()

#Set Monitors
pars1 <- c("alpha.hunt",
           "beta.econ.hunt",
           "H")

pars2 <- c("mu.econ.pred"
)

# Parallel Processing Setup
rm(out)
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
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

out <- clusterEvalQ(cl, {
  require(nimble)
  require(coda)
  model_test <- nimbleModel( code = code,
                             constants = constants,
                             data =  data,
                             inits = inits )
  
  model_test$simulate(c(
    'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'mu.econ.pred',
    'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt'
  ))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors = pars1, monitors2 =  pars2) 
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
save(files, file = 'model_output_HuntEff_pred.rdata')

