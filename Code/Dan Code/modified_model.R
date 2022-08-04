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
  
  ### Predator Predictor ###
  for(i in 1:4){
    sig.bbs.pred[i] ~ dgamma(1,1)
    mu.bbs.p[i] ~ dnorm(0, 1)
  }
  zeta.bbs[1] <- 1
  for(i in 2:4){
    zeta.bbs[i] ~ T(dnorm(0, 1),0,)
  }
  for(t in 1:n.year){
    pred.bbs.prime[t] ~ dnorm(0, 1) # Latent Predator Index
    ravens[t] ~ dnorm(mu.bbs.pred[1,t], sd = sig.bbs.pred[1])
    rthawk[t] ~ dnorm(mu.bbs.pred[2,t], sd = sig.bbs.pred[2])
    nharr[t] ~ dnorm(mu.bbs.pred[3,t], sd = sig.bbs.pred[3])
    pfal[t] ~ dnorm(mu.bbs.pred[4,t], sd = sig.bbs.pred[4])
    
    mu.bbs.pred[1,t] <- mu.bbs.p[1] + zeta.bbs[1] * pred.bbs.prime[t]
    mu.bbs.pred[2,t] <- mu.bbs.p[2] + zeta.bbs[2] * pred.bbs.prime[t]
    mu.bbs.pred[3,t] <- mu.bbs.p[3] + zeta.bbs[3] * pred.bbs.prime[t]
    mu.bbs.pred[4,t] <- mu.bbs.p[4] + zeta.bbs[4] * pred.bbs.prime[t]
  }
  
  #Winter Severity
  beta.awssi ~ dnorm(0, sd = 10)
  alpha.awssi ~ dnorm(0, sd = 10)
  for (r in 1:n.region) {
    sig.awssi[r] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for (t in 1:(n.year+1)) {
      awssi[r,t] ~ dnorm(alpha.awssi + (beta.awssi * era.awssi[t]), sd = sig.awssi[r])
    }
  }
  
  # Drought Index
  for(r in 1:n.region){
    sig.drought[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(t in 1:n.year){
      pdsi[t,r] ~ dnorm(0, sd = sig.drought[r])
    } #t
  } #r
  
  
  ################################################################################
  ### Hunter Effort ###
  mu.econ ~ dnorm(0, sd = 100)
  sig.econ ~ T(dt(0, pow(2.5, -2), 1), 0, )
  for (s in 1:n.species) {
    beta.econ.hunt[s] ~ dnorm(mu.econ, sd = sig.econ)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      sig.H[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      alpha.hunt[s,r] ~ dnorm(0, sd = 5) #sd = 1)
      for(k in 1:K){
        beta.spl.hunt[s,r,k] ~ dnorm(0, sd = sig.spl.hunt[s,r])
      } #k
      sig.spl.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      for(t in 1:n.year){ 
        latent.trend[s,t,r] <- inprod(beta.spl.hunt[s,r,1:K], basis[t,1:K])
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,t,r] <- alpha.hunt[s,r] + #intercept
          beta.econ.hunt[s] * pred.econ.prime[t] + #SEM economic indicator
          latent.trend[s,t,r]  #spline smoothing
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        H[s,t,r] <- exp(hunt.eps[s,t,r]) #Log Link
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        n.hunt[s,t,r] ~ dnorm(H[s,t,r], sd = sig.H[s,r]) #Number of hunters follows normal
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.hunt[s,t,r] <- H[s,t+1,r]/(H[s,t,r]+0.0001)
        log.r.hunt[s,t,r] <- log(lambda.hunt[s,t,r])
      } #t
    } #s
    
    for(t in 1:n.year){
      hunt.eps[1:n.species,t,r] ~ dmnorm(mu.hunt[1:n.species,t,r], cov =  Sigma.hunt[1:n.species,1:n.species,r] )
    } #t
    
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
  mu.bbs ~ dnorm(0, 0.01)
  sig.bbs ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.pdsi ~ dnorm(0, 0.01)
  sig.pdsi ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(r in 1:n.region){
    mu.hunter[r] ~ dnorm(0, 0.01)
    sig.hunter[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
  }
  
  for(s in 1:n.species){
    beta.wintsev.harv[s] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
    beta.bbs.harv[s] ~ dnorm(mu.bbs,  sd  = sig.bbs)
    beta.pdsi.harv[s] ~ dnorm(mu.pdsi, sd  = sig.pdsi)
  } #s
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.harv[s,r] ~ dnorm(0, sd = 5)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      beta.hunter.harv[s,r] ~ dnorm( mu.hunter[r], sd  =  sig.hunter[r])
      sig.N[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      # Process Model
      for(t in 1:n.year){
        #Unlinked estimate of Hunter Numbers
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        mu.harv[s,t,r] <- alpha.harv[s,r] +                       # intercepts
          beta.hunter.harv[s,r] * latent.trend[s,t,r] + # Latent Hunter Trend
          beta.wintsev.harv[s] * awssi[r,t]  +         # Previous winter severity (Affecting Survival)
          beta.pdsi.harv[s] * pdsi[t,r] +           # Same year, spring/summer drought (Affecting Survival/Reproduction)
          beta.bbs.harv[s] * pred.bbs.prime[t]     # Latent predator index (Affecting Reproduction)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        N[s,t,r] <- exp(harv.eps[s,t,r]) #Log Link
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        n.harv[s,t,r] ~ dnorm(N[s,t,r], sd = sig.N[s,r]) #Number of hunters follows normal
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
  ### Chukar Site Abundance ###
  for(r in 1:n.region){
    theta.chuk[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    mod.chuk[r] ~ dlogis(0,1)
    mu.chuk[r] ~ dlogis(0,1)
  }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.chuk){
      # Change the latent trend number (3) to match chukar
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      log.r.chuk[p,t-1] <-  mu.chuk[reg.chuk[p]] + mod.chuk[reg.chuk[p]] * latent.trend[3, t+13, reg.chuk[p]] #log.r.harv[t=13] is 1990-1991
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
      C.chuk[p,t] <- exp(log.r.chuk[p,t-1]) * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[reg.chuk[p]]/(theta.chuk[reg.chuk[p]] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(prob = rate.chuk[p,t-1], size = theta.chuk[reg.chuk[p]]) #obs. # of chukars follow neg-bin
    } #t
  } #p
  
  ################################################################################
  ### Birds per Hunter
  for(t in 1:n.year){
    for(s in 1:n.species){
      for(r in 1:n.region){
        BPH[s,t,r] <-  N[s,t,r]/H[s,t,r]
        BPH2[s,t,r]<- exp(mu.hunt[s,t,r]-mu.harv[s,t,r])
      }
    }
  }
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
  une = une, #BL Unemployment information for Nevada, scaled
  res = res, #Residential license sales
  pdi = scale(PDI)[,1], #personal disposable income
  gas = scale(GAS)[,1], #gas prices
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  # wpdsi = wpdsi, #winter drought index, scaled
  # rabbits = rabbits, #Number of rabbits harvested
  # raven = as.vector(bbs.df$raven)[-nrow(bbs.df)], #bbs bayes index for ravens, t = 1 is 1975
  # nharrier = as.vector(bbs.df$nharrier), #bbs bayes index for northern harriers, t = 1 is 1975
  
  ### Hunter Effort
  n.hunt = hunters/1000, #Observed number of hunters for each species each year
  basis = B,
  ravens = bbs.df[,2], rthawk = bbs.df[,3], nharr =  bbs.df[,4], pfal = bbs.df[,5],
  ### Total Harvest
  n.harv = upland/1000,
  
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
  K = dim(B)[2],
  
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
Ni[,1,] <- (upland[,1,] + 50)/1000

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
#C.chuk.init[,2:ncol(C.chuk.init)] <- NA

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
  
  mu.chuk = c(0,0),
  ### Predictors
  #Winter Severity
  sig.awssi = rep(1,2),
  beta.awssi = 0,
  alpha.awssi = 0,
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
  mu.hunter = c(0,0),
  sig.hunter = c(1,1),
  beta.wintsev.harv = rep(0, 7),
  beta.pdsi.harv = rep(0, 7),
  beta.bbs.harv = rep(0, 7),
  beta.hunter.harv = matrix(0, 7,2),
  beta.econ.hunt = rep(0, 7),
  beta.spl.hunt = array(0, dim = c(7,2,dim(B)[2])),
  sig.N =  matrix(1, 7,2),
  mu.wintsev.hunt = 0,
  sig.wintsev.hunt = 1,
  mu.jobs = 0,
  sig.jobs = 1,
  beta.jobs = rep(0, 7),
  beta.wintsev.hunt = rep(0, 7),
  
  ### Hunter Effort
  sig.H =  matrix(1, 7,2),
  n.hunt = n.hunt.i/1000,
  Q.hunt = abind(Q,Q,along = 3),
  P.hunt = abind(P,P,along = 3),
  Lambda.hunt = abind(diag(n.species),diag(n.species),along = 3),
  Delta.hunt = abind(Delta,Delta,along = 3),
  rho.hunt = abind(diag(n.species),diag(n.species),along = 3),
  alpha.hunt = matrix(0, ncol = 2, nrow = 7),
  sig.hunt = matrix(1, ncol = 2, nrow = 7),
  sig.spl.hunt = matrix(1, ncol = 2, nrow = 7),
  sig.econ.pred = c(1,1,1,1),
  ### Total Harvest
  n.harv = n.harv.i/1000,
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
  n.chuk = as.matrix(C.chuk.init),
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
  #"beta.wintsev.hunt",
  #"beta.jobs",
  #'zeta','pred.prime',
  'beta.spl.hunt',
  ### Total Harvest
  "alpha.harv",'sig.harv','harv.eps',
  'beta.pdsi.harv',
  "beta.hunter.harv"
  #"beta.raven.harv"
)

pars2 <- c(### Hunter Effort
  "H",'sig.H',
  "rho.hunt",'rho.harv',
  "log.r.hunt",
  ### Total Harvest
  "N",'sig.N',
  "log.r.harv",
  # ### Chukar Site Abundance
  "log.r.chuk",'mod.chuk',
  ### Birds per Hunter
  "BPH",'BPH2')

pars3 <- c(
  'pred.econ.prime','pred.bbs.prime','sig.econ.pred', 'sig.bbs.pred',
  'mu.bbs.p','zeta.econ','zeta.bbs','alpha.awssi', 'beta.awssi','awssi',
  'beta.econ.hunt','sig.H','alpha.hunt','beta.spl.hunt','latent.trend',
  'beta.wintsev.harv','beta.bbs.harv','beta.pdsi.harv','beta.hunter.harv','sig.N','mu.harv','mu.hunt','mod.chuk'
)

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

clusterExport(cl, c("code", "inits", "data", "constants", "pars1", "pars2",'pars3', "pars.pred")) #identify what is to be exported to each cluster

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
    'mu.econ.p', 'zeta.econ', 'pred.econ.prime', 'gas', 'res', 'pdi', 'mu.econ.pred', 'sig.bbs.pred', 'mu.bbs.p', 'zeta.bbs',
    'pred.bbs.prime', 'mu.bbs.pred', 'sig.drought', 'mu.econ', 'sig.econ', 'latent.trend', 
    'mu.hunt', 'H', 'n.hunt', 'lambda.hunt', 'log.r.hunt', 'hunt.eps', 'mu.bbs', 'sig.bbs', 'mu.pdsi', 'sig.pdsi',
    'mu.harv', 'N', 'lambda.harv', 'log.r.harv', 'harv.eps', 'C.chuk', 'log.r.chuk', 'rate.chuk', 'BPH','BPH2'
  ))
  
  model_test$initializeInfo()
  model_test$calculate()
  
  mcmcConf <-  configureMCMC( model_test,   monitors =  c(pars1, pars2,pars3))
  # mcmcConf <-  configureMCMC( model_test,   monitors =  pars1, monitors2 = pars2)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 75000, niter = 125000, thin = 5, thin2 = 5)
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


beta.hunter.harv    <- MCMCsummary(mcmcList1, params = 'beta.hunter.harv')
beta.wintsev.harv   <- MCMCsummary(mcmcList1, params = 'beta.wintsev.harv')
beta.pdsi.harv      <- MCMCsummary(mcmcList1, params = 'beta.pdsi.harv')
beta.bbs.harv       <- MCMCsummary(mcmcList1, params = 'beta.bbs.harv')
latent.trend        <- apply(MCMCpstr(mcmcList1, params = 'latent.trend',type='chains')$latent.trend,c(1:3), mean)
latent.trend_se     <- apply(MCMCpstr(mcmcList1, params = 'latent.trend',type='chains')$latent.trend,c(1:3), sd)
alpha.hunt          <- MCMCsummary(mcmcList1, params = 'alpha.hunt')
beta.econ.hunt      <- MCMCsummary(mcmcList1, params = 'beta.econ.hunt')

n <- 1000
bph_prediction <- hunt_prediction <- harv_prediction <- array(NA, dim = c(7,2,n))

cov <- rnorm(4, 0, 1)

for(i in 1:7){
  harv_prediction[i,1,] <- rnorm(n, alpha.harv$mean[i], alpha.harv$sd[i]) + 
    rnorm(n, beta.hunter.harv$mean[i], beta.hunter.harv$sd[i]) * rnorm(n, latent.trend[i,42,1], latent.trend_se[i,42,1]) +
    rnorm(n, beta.wintsev.harv$mean[i], beta.wintsev.harv$sd[i]) * cov[1] +
    rnorm(n, beta.pdsi.harv$mean[i], beta.pdsi.harv$sd[i]) * cov[2] +
    rnorm(n, beta.bbs.harv$mean[i], beta.bbs.harv$sd[i]) * cov[3]
  
  harv_prediction[i,2,] <- rnorm(n, alpha.harv$mean[i+7], alpha.harv$sd[i+7]) + 
    rnorm(n, beta.hunter.harv$mean[i+7], beta.hunter.harv$sd[i+7]) * rnorm(n, latent.trend[i,42,2], latent.trend_se[i,42,2]) +
    rnorm(n, beta.wintsev.harv$mean[i], beta.wintsev.harv$sd[i])* cov[1]+
    rnorm(n, beta.pdsi.harv$mean[i], beta.pdsi.harv$sd[i]) * cov[2] +
    rnorm(n, beta.bbs.harv$mean[i], beta.bbs.harv$sd[i]) * cov[3]
  
  
  hunt_prediction[i,1,] <-rnorm(n, alpha.hunt$mean[i], alpha.hunt$sd[i]) + 
    rnorm(n, latent.trend[i,42,1], latent.trend_se[i,42,1]) +
    rnorm(n, beta.econ.hunt$mean[i], beta.econ.hunt$sd[i]) * cov[4]
  
  hunt_prediction[i,2,] <-rnorm(n, alpha.hunt$mean[i+7], alpha.hunt$sd[i+7]) + 
    rnorm(n, latent.trend[i,42,2], latent.trend_se[i,42,2]) +
    rnorm(n, beta.econ.hunt$mean[i], beta.econ.hunt$sd[i]) * cov[4]
}


bph_prediction <- exp(harv_prediction - hunt_prediction)
# Number of harvest birds (in thousands)
exp(apply(harv_prediction, c(1:2), quantile, probs = c(0.075, .5, 0.935) ))
# Number of huntres (in thousands)
exp(apply(hunt_prediction, c(1:2), quantile, probs = c(0.075, .5, 0.935)  ))
# Birds per hunter
apply(bph_prediction, c(1:2), quantile,probs = c(0.075, .5, 0.935)  )







save(mcmcList1, file = 'modifed_model.rdata')

test.bph2 <- MCMCsummary(mcmcList1, 'BPH2', hpd_prob = 0.85, HPD = TRUE)
test.bph <- MCMCsummary(mcmcList1, 'BPH', hpd_prob = 0.85, HPD = TRUE)
test.H   <- MCMCsummary(mcmcList1, 'H', hpd_prob = 0.85, HPD = TRUE)
test.N  <- MCMCsummary(mcmcList1, 'N', hpd_prob = 0.85, HPD = TRUE)


hunt_2017 <- subset(harvest_full, Year == 2017)
hunters_east <- c(377,55,2305,329,25,116,820)
hunters_west <- c(241,1629,8352,210,61,394,500)

h.prime <- hunters
h.prime[1:7,42,1] <- hunters_east
h.prime[1:7,42,2] <- hunters_west

animals_east <- c(605,225,17124,1485,78,419,1522)
animals_west <- c(320,11617,57391,646,98,1757,812)

u.prime <- upland
u.prime[1:7,42,1] <- animals_east
u.prime[1:7,42,2] <- animals_west

test.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                        year = 1976:2017, 
                                        region = c('west','east') ),
                            test.H, hunters = c(h.prime)/1000)

test.dfn <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                         year = 1976:2017, 
                                         region = c('west','east') ),
                             test.N, upland = c(u.prime)/1000)

test.dfbp <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                          year = 1976:2017, 
                                          region = c('west','east') ),
                              test.bph, upland = c(u.prime), hunters = c(h.prime))

test.dfbp2 <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                           year = 1976:2017, 
                                           region = c('west','east') ),
                               test.bph2, upland = c(u.prime), hunters = c(h.prime))


fig_upland <-   ggplot(data = test.dfn, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `85%_HPDL`,ymax = `85%_HPDU`), linetype = 'dashed') +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = (1), color = NA) +  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested', x = 'Year') +
  scale_fill_flat_d() + theme_modern() +
  geom_point(aes(x = year, y = (upland)), color = 'red') 
fig_upland

fig_succ <-   ggplot(data = test.dfbp, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `85%_HPDL`,ymax = `85%_HPDU`),   linetype = 'dashed') +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested per hunter', x = 'Year') +
  scale_fill_flat_d() + theme_modern() +
  geom_point(aes(x = year, y = upland/hunters), color = 'red') 
fig_succ

fig_succ1 <-   ggplot(data = subset(test.dfbp, year == 2017), aes(x = species, y = mean)) +
  geom_pointrange2(aes(x = species, y = mean, ymin = `85%_HPDL`,ymax = `85%_HPDU`,fill = region), size = .66, shape = 21, linetype = 'dashed', position = position_dodge2(.5)) +
  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested per hunter\nduring missing year', x = 'Species', fill = 'Region') +
  scale_fill_flat_d() + theme_modern() +
  geom_point(aes(x = species, y = upland/hunters), fill = 'gray',shape = 21, size = 3, position = position_dodge2(.5)) 
fig_succ1


ggplot(data = test.df, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean,ymin = `85%_HPDL`,ymax = `85%_HPDU`), fill = 'black', shape = 21, size = .75, position = position_dodge(1)) +
  geom_point(aes(x = year, y = (hunters)), color = 'red') +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1) +
  scale_fill_flat_d() + theme_modern()





