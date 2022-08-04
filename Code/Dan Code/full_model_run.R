### Modeled Variation in N instead of log.r.harv
#### Spline to deal with serial autocorrelation in N/H
#### No Covariates

lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2016 #Last year from which data will be used
final.y <- 2017 #Last year to predict
year.hold <- cutoff.y +1

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- final.y - cutoff.y
load("./model_output_HuntEff_pred.rdata")


### Model Code
code <- nimbleCode( {
  # ################################################################################
  # ### Predictors ###
  
  for(i in 1:4){
    sig.pred[i] ~ dgamma(1,1)
        mu_p[i] ~ dnorm(0, 1)
  }
  zeta[1] <- 1
  for(i in 2:4){
    zeta[i] ~ dnorm(0, 1)
  }
  for(t in 1:n.year){
    pred.prime[t] ~ dnorm(0, 1) # Latent Predator Index
        ravens[t] ~ dnorm(mu.pred[1,t], sd = sig.pred[1])
        rthawk[t] ~ dnorm(mu.pred[2,t], sd = sig.pred[2])
         nharr[t] ~ dnorm(mu.pred[3,t], sd = sig.pred[3])
          pfal[t] ~ dnorm(mu.pred[4,t], sd = sig.pred[4])
    
      mu.pred[1,t] <- mu_p[1] + zeta[1] * pred.prime[t]
      mu.pred[2,t] <- mu_p[2] + zeta[2] * pred.prime[t]
      mu.pred[3,t] <- mu_p[3] + zeta[3] * pred.prime[t]
      mu.pred[4,t] <- mu_p[4] + zeta[4] * pred.prime[t]
  }


  ################################################################################
  ### Hunter Effort ###
  mu.wintsev.hunt ~ dnorm(0, sd = 100)
  sig.wintsev.hunt ~ T(dt(0, pow(2.5, -2), 1), 0, )
  mu.jobs ~ dnorm(0, sd = 100)
  sig.jobs ~ T(dt(0, pow(2.5, -2), 1), 0, )
  
  for (s in 1:n.species) {
    beta.wintsev.hunt[s] ~ dnorm(mu.wintsev.hunt, sd = sig.wintsev.hunt)
    beta.jobs[s] ~ dnorm(mu.jobs, sd = sig.jobs)
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
          #beta.wintsev.hunt[s] * awssi[r, t + 1] + 
          beta.jobs[s] * une[t] + 
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
  
  
  ################################################################################
  ### Total Harvest ###
  mu.wintsev.harv ~ dnorm(0, 0.01)
  sig.wintsev.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.pdsi.harv ~ dnorm(0, 0.01)
  sig.pdsi.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.raven ~ dnorm(0, 0.01)
  sig.raven ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.hunter ~ dnorm(0, 0.01)
  sig.hunter ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
    beta.wintsev.harv[s] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
       beta.pdsi.harv[s] ~ dnorm(mu.pdsi.harv, sd  = sig.pdsi.harv)
      beta.raven.harv[s] ~ dnorm(mu.raven, sd  = sig.raven)
     beta.hunter.harv[s] ~ dnorm(mu.hunter, sd  = sig.hunter)
  } #s

  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.harv[s,r] ~ dnorm(5, sd = 3)
      for(k in 1:K){
        beta.spl.harv[s,r,k] ~ dnorm(0, sd = sig.spl.harv[s,r])
      } #k
      sig.spl.harv[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      for(t in 1:n.year){
        #Unlinked estimate of Hunter Numbers
        mu.harv[s,t,r] <- alpha.harv[s,r] +#regression formula
                          beta.hunter.harv[s] * ((n.hunt[s,t,r] - mean(n.hunt[s,1:n.year,r]))/sd(n.hunt[s,1:n.year,r])) + #Current season hunter numbers
                          beta.wintsev.harv[s] * awssi[r,t]  + #Previous winter severity (Affecting Survival)
                          beta.pdsi.harv[s] * pdsi[t,r] +
                          beta.raven.harv[s] * pred.prime[t] + #Previous Year BBS index (Affecting Reproduction)
                          inprod(beta.spl.harv[s,r,1:K], Z.harv[t,1:K,s,r]) #spline smoothing
        
        N[s,t,r] <- exp(harv.eps[s,t,r]) #Log Link
        
        n.harv[s,t,r] ~ dpois(N[s,t,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.harv[s,t,r] <- N[s,t+1,r]/N[s,t,r]
        log.r.harv[s,t,r] <- log(lambda.harv[s,t,r])
      }
    } #s

    for(t in 1:n.year){
      harv.eps[1:n.species,t,r] ~ dmnorm(mu.harv[1:n.species,t,r], cov =  Sigma.harv[1:n.species,1:n.species,r] )
    }
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
  
  ################################################################################
  ### Birds per Hunter
  for(t in 1:n.year){
    for(s in 1:n.species){
      for(r in 1:n.region){
        BPH[s,t,r] <- N[s,t,r]/H[s,t,r]
      }
    }
  }
  
  ################################################################################
  ### Chukar Site Abundance ###
  for(r in 1:n.region){
    theta.chuk[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    mod.chuk[r] ~ dlogis(0,1)
  }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.chuk){
      log.r.chuk[p,t-1] <-  mod.chuk[reg.chuk[p]] * log.r.harv[3, t+13, reg.chuk[p]] #log.r.harv[t=13] is 1990-1991
      
      C.chuk[p,t] <- exp(log.r.chuk[p,t-1]) * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[reg.chuk[p]]/(theta.chuk[reg.chuk[p]] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(prob = rate.chuk[p,t-1], size = theta.chuk[reg.chuk[p]]) #obs. # of chukars follow neg-bin
    } #t
  } #p
  
})


### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2016 #Last year from which data will be used
final.y <- 2017 #Last year to predict
year.hold <- cutoff.y +1

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- final.y - cutoff.y
source("./ChukSEM - Data Prep.R")


### Run Hunter Effort Solo Model to get estimates of H to create spline inputs
# Sys.time()
# source("./ChukSEM - Hunter Effort Model - Predict.R")
# Sys.time()

load("./model_output_HuntEff_pred.rdata")
# mcmcList2 <- files[[1]]

### Load Model Code
source("./ChukSEM - Model Only - Base w Covs.R")


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
  ### Covariates
  une = une, #BL Unemployment information for Nevada, scaled
  # res = res, #Residential license sales
  # PDI = PDI, #personal disposable income
  # GAS = GAS, #gas prices
  awssi = awssi, #winter severity index, scaled
    pdsi = pdsi, #Previous breeding season drought index
  # wpdsi = wpdsi, #winter drought index, scaled
  # rabbits = rabbits, #Number of rabbits harvested
  # raven = as.vector(bbs.df$raven)[-nrow(bbs.df)], #bbs bayes index for ravens, t = 1 is 1975
  # nharrier = as.vector(bbs.df$nharrier), #bbs bayes index for northern harriers, t = 1 is 1975
  
  ### Hunter Effort
  n.hunt = hunters, #Observed number of hunters for each species each year
  Z.hunt = ZZ1, #Spline
  ravens = bbs.df[,2], rthawk = bbs.df[,3], nharr =  bbs.df[,4], pfal = bbs.df[,5],
  ### Total Harvest
  n.harv = upland,
  Z.harv = ZZ, #Spline
  
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
  era.awssi = c(rep(0,length(1975:1994)),rep(1, length(1995:2001)), rep(0, length(2002:2017))), #Groupings for change in gas prices
  
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
  
  zeta = c(1,1,1,1), sig.pred = c(1,1,1,1),mu_p = c(0,0,0,0), pred.prime = rnorm(ncol(hunters)),zeta.upland = matrix(0, 7,2),
  
  
  ### Predictors
  #Winter Severity
  sig.awssi = rep(1,2),
  beta.awssi = 0,
  alpha.awssi = 0,
  awssi = awssi.init,
  # #Ravens
  alpha.rav = 0,
  beta.t.rav = 0,
  sig.rav = 1,
  ar1.rav = 0,
  #Unemployment
  sig.une = 1,
  une = une.init,
  
  ###Beta Coefficient Means/SDs
  mu.wintsev.harv = 0,
  sig.wintsev.harv = 1,
  mu.pdsi.harv = 0,
  sig.pdsi.harv = 1,
  mu.raven = 0,
  sig.raven = 1,
  mu.hunter = 0,
  sig.hunter = 1,
  beta.wintsev.harv = rep(0, 7),
  beta.pdsi.harv = rep(0, 7),
  beta.raven.harv = rep(0, 7),
  beta.hunter.harv = rep(0, 7),
  beta.spl.harv = array(0, dim = c(7,2,12)),
  mu.wintsev.hunt = 0,
  sig.wintsev.hunt = 1,
  mu.jobs = 0,
  sig.jobs = 1,
  beta.jobs = rep(0, 7),
  beta.wintsev.hunt = rep(0, 7),
  
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
  'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt','zeta','pred.prime',
  'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',
  'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
  'BPH'
))
model_test$initializeInfo()
model_test$calculate()


pars1 <- c(### Hunter Effort
  "alpha.hunt",
  "beta.wintsev.hunt",
  "beta.jobs",
  'zeta','pred.prime',
  'beta.spl.hunt',
  ### Total Harvest
  "alpha.harv",'sig.harv','harv.eps',
  'beta.pdsi.harv',
  'beta.spl.harv',
  "beta.wintsev.harv",
  "beta.hunter.harv",
  "beta.raven.harv"
)

pars2 <- c(### Hunter Effort
  "H",
  "rho.hunt",'rho.harv',
  "log.r.hunt",
  ### Total Harvest
  "N",
  "log.r.harv",
  # ### Chukar Site Abundance
  "log.r.chuk",
  ### Birds per Hunter
  "BPH")

pars.pred <- c(
  "sig.une",
  "une",
  "alpha.rav",
  "beta.t.rav",
  "sig.rav",
  "ar1.rav",
  "raven",
  "sig.awssi",
  "beta.awssi",
  "alpha.awssi",
  "awssi"
)

# Parallel Processing Setup
rm(out)
rm(out.full.predict)
start_time <- Sys.time() # To track runtime
start_time
nc <- 3   # number of chains
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
    'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt', 'log.r.hunt','zeta','pred.prime',
    'beta.spl.harv', 'pred.spl.harv','mu.harv', 'harv.eps', 'N', 'Sigma.harv', 'lambda.harv', 'log.r.harv',
    'theta.chuk','rate.chuk', 'log.r.chuk', 'C.chuk', 'mod.chuk', 'chuk.eps',
    'BPH'
  ))
  
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors =  c(pars1, pars2))
  # mcmcConf <-  configureMCMC( model_test,   monitors =  pars1, monitors2 = pars2)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)

  samplesList <- runMCMC(Cmcmc,nburnin = 75000, niter = 150000, thin = 5, thin2 = 5)
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

#Save Outputs as file
files <- list(mcmcList1, code)
# files <- list(mcmcList1, mcmcList2, code)
save(files, file = 'hold.rdata')


test.bph <- MCMCsummary(mcmcList1, 'BPH')
test.H   <- MCMCsummary(mcmcList1, 'H')
test.N   <- MCMCsummary(mcmcList1, 'N')


hunt_2017 <- subset(harvest_full, Year == 2017)


hunters_east <- c(377,55,2305,329,25,116,820)
hunters_west <- c(241,1629,8352,210,61,394,500)

hunters[1:7,42,1] <- hunters_east
hunters[1:7,42,2] <- hunters_west

animals_east <- c(605,225,17124,1485,78,419,1522)
animals_west <- c(320,11617,57391,646,98,1757,812)

upland[1:7,42,1] <- animals_east
upland[1:7,42,2] <- animals_west

test.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                        year = 1976:2017, 
                                        region = c('west','east') ),
                            test.H, hunters = c(hunters))

test.dfn <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                         year = 1976:2017, 
                                         region = c('west','east') ),
                             test.N, upland = c(upland))


test.dfbp <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                         year = 1976:2017, 
                                         region = c('west','east') ),
                              test.bph, upland = c(upland), hunters = c(hunters))



fig_upland <-   ggplot(data = test.dfn, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`),  size = .75) +
  geom_pointrange2(data = subset(test.dfn,year>2017), aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), color = 'red',  size = .75) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1000, color = NA) +  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested', x = 'Year') +
  scale_fill_flat_d() + theme_modern() +
  geom_point(aes(x = year, y = upland), color = 'red') 
fig_upland

fig_succ <-   ggplot(data = test.dfbp, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`),  size = .75) +
  geom_pointrange2(data = subset(test.dfbp,year>2017), aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), color = 'red',  size = .75) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested per hunter', x = 'Year') +
  scale_fill_flat_d() + theme_modern() +
  geom_point(aes(x = year, y = upland/hunters), color = 'red') 
fig_succ

ggplot(data = test.df, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), fill = 'black', shape = 21, size = .75, position = position_dodge(1)) +
  geom_point(aes(x = year, y = hunters), color = 'red') +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1) +
  scale_fill_flat_d() + theme_modern()

