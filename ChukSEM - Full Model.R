### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Run Hunter Effort Solo Model to get estimates of H to create spline inputs
source("./ChukSEM - Hunter Effort Model.R")

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Predictors
  #Personal Disposable Income
  alpha.pdi ~ dnorm(0, 0.001) #intercept
  beta.t.pdi ~ dnorm(0, 0.01) #year
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1 ~ dunif(-1,1) #Autoregressive parameter
  
  pdi.trend[1] <- alpha.pdi + beta.t.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    pdi.trend[t] <- alpha.pdi + beta.t.pdi * t
    mu.pdi[t] <- pdi.trend[t] + ar1 * (PDI[t-1] - pdi.trend[t-1])
  } #t
  
  #Gas Prices
  alpha.gas ~ dnorm(1.5, 1)
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.gas[1] ~ dnorm(0, 0.01)
  beta.gas[2] ~ dnorm(0, 0.01)
  
  #Relative Cost of Gas
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(alpha.gas + beta.gas[era[t]]*t, sd = sig.gas),0,)
    rel.cost[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  } #t
  
  #Unemployment Rate
  sig.une ~ dunif(0,5)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
  } #t
  
  #Drought Index
  for(r in 1:n.region){
    sig.wpdsi[r] ~ dunif(0,5)
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
    } #t
  } #r
  
  ################################################################################
  ### Hunter Effort ###
  for(s in 1:n.species){
    beta.drought.hunt[s] ~ dnorm(0, sd = 100)
    beta.wintsev.hunt[s] ~ dnorm(0, sd  = 100)
    beta.jobs[s] ~ dnorm(0, sd = 100)
    beta.income[s] ~ dnorm(0, sd  = 100)
    beta.license[s] ~ dnorm(0, sd  = 100)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.hunt[s,r] ~ dnorm(5, sd = 1)
      for(k in 1:K){
        beta.spl.hunt[s,r,k] ~ dnorm(0, sd = sig.spl.hunt[s,r])
      } #k
      sig.spl.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(t in 1:n.year){ 
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,r,t] <- alpha.hunt[s,r] + #intercept
          beta.drought.hunt[s] * wpdsi[t,r] + #concurrent winter drought index
          beta.wintsev.hunt[s] * awssi.hunt[r,t] + #concurrent winter severity
          beta.jobs[s] * une[t] + #concurrent years unemployment
          beta.income[s] * rel.cost[t] + #PDI/Gas Price
          beta.license[s] * res[t] + #Hunting licences sold that season
          inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #spline smoothing
        
        pred.spl.hunt[s,r,t] <- inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #Derive spline smoothing for examination later
        
        log(H[s,t,r]) <- hunt.eps[s,r,t] #Log Link
        n.hunt[s,t,r] ~  dpois(H[s,t,r]) #Number of hunters follows Poisson
      } #t
    } #s
    
    for(t in 1:n.year){
      hunt.eps[1:n.species,r,t] ~ dmnorm(mu.hunt[1:n.species,r,t], cov =  Sigma.hunt[1:n.species,1:n.species,r] )
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
  for(s in 1:n.species){
    beta.drought.harv[s] ~ dnorm(0, sd = 100)
    beta.wintsev.harv[s] ~ dnorm(0, sd  = 100)
    beta.rabbit[s] ~ dnorm(0, sd  = 100)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      # Regression Coefficients
      alpha.harv[s,r] ~ dnorm(0, sd = 1)
      for(k in 1:K){
        beta.spl.harv[s,r,k] ~ dnorm(0, sd = sig.spl.harv[s,r])
      } #k
      sig.spl.harv[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Total harvest, Year 1
      
      for(t in 2:(n.year)){
        mu.harv[s,t-1,r] <- alpha.harv[s,r] +#regression formula
          beta.drought.harv[s] * pdsi[t,r] + #previous breeding season drought index
          beta.wintsev.harv[s] * awssi[r,t] + #concurrent winter severity
          beta.rabbit[s] * rabbits[t,r] + #concurrent rabbit harvest
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
  
  ################################################################################
  ### Sage Grouse Wing-Bee ###
  ### Chukar Site Abundance
  beta.drought.sg ~ dnorm(0, 0.01)
  beta.wintsev.sg ~ dnorm(0, 0.01)
  beta.rabbit.sg ~ dnorm(0, 0.01)
  
  for(r in 1:n.region){
    alpha.sg[r] ~ dnorm(0, sd = 100) #Intercept for
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.sg[r,1] ~ dpois(wing.b[r,1]) #Equivalent of Poisson lambda
    
    sigma.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod.sg[r] ~ dlogis(0,1)
    
    for(t in 2:n.years.sg){
      mu.sg[r, t-1] <- mod.sg[r] * log.r.harv[7, t+24, r] #log.r.harv[t=29] is 2004-2005
      sg.eps[r,t-1] ~ dnorm(mu.sg[r, t-1], sd = sigma.sg[r])
      
      log.r.sg[r,t-1] <- alpha.sg[r] + 
        sg.eps[r,t-1] + #Variation associated with change in harvest
        beta.drought.sg * pdsi[t+25,r] + #previous breeding season drought index
        beta.wintsev.sg * awssi[r,t+25] + #previous year's winter severity
        beta.rabbit.sg * rabbits[t+25,r] #concurrent rabbit harvest
      
      lambda.sg[r,t-1] <- exp(log.r.sg[r,t-1]) # Tranform between finite and instantaneous growth rate
      
      C.sg[r,t] <- lambda.sg[r,t-1] * C.sg[r,t-1] #Change in Count over time
      
      rate.sg[r,t-1] <- theta.sg[r]/(theta.sg[r] + C.sg[r,t]) #NB rate
      wing.b[r,t] ~ dnegbin(rate.sg[r,t-1], theta.sg[r]) #Wing-Bee data follows negative binomial
    } #t
  } #r
  
  ################################################################################
  ### Chukar Site Abundance ###
  beta.drought.chuk ~ dnorm(0, 0.01)
  beta.wintsev.chuk ~ dnorm(0, 0.01)
  beta.rabbit.chuk ~ dnorm(0, 0.01)
  
  for(p in 1:n.site){
    alpha.chuk[p] ~ dnorm(0, sd = 100) #Intercept
    
    theta.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod.chuk[p] ~ dlogis(0,1)
    
    for(t in 2:n.year.chuk){
      mu.chuk[p, t-1] <- mod.chuk[p] * log.r.harv[3, t+13, chuk.reg[p]] #log.r.harv[t=14] is 1990-1991
      chuk.eps[p,t-1]  ~ dnorm(mu.chuk[p, t-1], sd = sigma.chuk[p]) #Coeffient for annual change
      
      log.r.chuk[p,t-1] <- alpha.chuk[p] + 
        chuk.eps[p,t-1] + #Unlinked change in abundance, log.harv.chuk[t=1] is 1990-1991
        beta.drought.chuk * pdsi[t+13,chuk.reg[p]] + #previous breeding season drought index
        beta.wintsev.chuk * awssi[chuk.reg[p],t+14] + #previous year's winter severity
        beta.rabbit.chuk * rabbits[t+13,chuk.reg[p]] #previous year's rabbit harvest
      lambda.chuk[p,t-1] <- exp(log.r.chuk[p,t-1]) #Change in abundance
      
      C.chuk[p,t] <- lambda.chuk[p,t-1] * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[p]/(theta.chuk[p] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(rate.chuk[p,t-1], theta.chuk[p]) #obs. # of chukars follow neg-bin
    } #t
  } #p 
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

# Time?
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
  res = scale(res)[,1], #Residential license sales
  PDI = PDI, #personal disposable income
  GAS = GAS, #gas prices
  awssi = awssi, #winter severity index, scaled
  awssi.hunt = awssi[,-1], #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  wpdsi = wpdsi, #winter drought index, scaled
  rabbits = rabbits, #Number of rabbits harvested 
  
  ### Hunter Effort
  n.hunt = hunters, #Observed number of hunters for each species each year
  Z.hunt = ZZ1, #Spline 
  
  ### Total Harvest
  n.harv = upland,
  Z.harv = ZZ, #Spline
  
  ### Sage Grouse WingBee
  wing.b =  wing.b,
  
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
  n.species = 7,
  n.year = ncol(hunters),
  K = 12,
  
  ### Predictors
  era = c(rep(1,19),rep(2, 27)), #Groupings for change in gas prices 
  
  ### Hunter Effort
  I.hunt = abind(I,I,along = 3),
  
  ### Total Harvest
  I.harv = abind(I2,I2,along = 3),
  
  ### Sage Grouse Wing-Bee
  n.years.sg = n.years.sg,
  
  ### Chukar Site Abundance
  n.site = 13,
  n.year.chuk = ncol(chukar),
  chuk.reg = chuk.reg
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

PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), NA)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, NA)

##Total Harvest
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
C.sg.i[,1] <- floor(rowMeans(wing.b, na.rm = T))

## Chukar Site Abundance
chukar_na <- chukar 
chukar_na <- ifelse(is.na(chukar == TRUE), floor(mean(as.matrix(chukar), na.rm = T)), NA)


# Wrapper Function
initsFunction <- function() list(
  ### Predictors
  GAS = GAS.inits,
  PDI = PDI.inits,
  alpha.pdi = 2.9,
  beta.t.pdi = 0,
  sig.pdi = 1,
  ar1 = 0,
  alpha.gas = .9,
  beta.gas = rep(0,2),
  sig.gas = 1,
  sig.une = 1,
  sig.wpdsi = rep(1,2),

  sig.spl.hunt = matrix(1, ncol = 2, nrow = 7),
  # mu.drought.hunt = rep(0,2),
  # sig.drought.hunt = rep(1,2),
  # mu.wintsev.hunt = rep(0,2),
  # sig.wintsev.hunt = rep(1,2),
  # mu.jobs = rep(0,2),
  # sig.jobs = rep(1,2),
  # mu.income = rep(0,2),
  # sig.income = rep(1,2),
  # mu.license = rep(0,2),
  # sig.license = rep(1,2),
  
  sig.spl.harv = matrix(1, ncol = 2, nrow = 7),
  # mu.drought.harv = rep(0,2),
  # sig.drought.harv = rep(1,2),
  # mu.wintsev.harv = rep(0,2),
  # sig.wintsev.harv = rep(1,2),
  # mu.rabbit = rep(0,2),
  # sig.rabbit = rep(1,2),
  
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
  beta.drought.hunt = rep(0, 7),
  beta.wintsev.hunt = rep(0, 7),
  beta.jobs = rep(0, 7),
  beta.income = rep(0, 7),
  beta.license = rep(0, 7),
  
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
  beta.drought.harv = rep(0, 7),
  beta.wintsev.harv = rep(0, 7),
  beta.rabbit = rep(0, 7),
  beta.spl.harv = array(0, dim = c(7,2,12)),
  
  ### Sage Grouse Wing-Bee
  theta.sg = rep(1,2),
  # sigma.sg = rep(0,2),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  C.sg = C.sg.i,
  # mod.sg = rep(1,2),
  alpha.sg =  rep(0,2),
  beta.drought.sg = 0,
  beta.wintsev.sg = 0,
  beta.rabbit.sg = 0,
  # alpha.sg =  rep(0,2),
  # beta.drought.sg = rep(0, 2),
  # beta.wintsev.sg = rep(0, 2),
  # beta.rabbit.sg = rep(0, 2),
  sigma.sg = rep(1,2),
  sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  mod.sg = rep(1,2),
  
  ### Chukar Site Abundance
  theta.chuk = rep(1,13),
  # sigma.chuk = rep(0,2),
  # sg.eps = matrix(0, nrow = 2, ncol = n.years.chuk-1),
  n.chuk = chukar_na,
  # mod.chuk = rep(1,2),
  alpha.chuk =  rep(0,13),
  beta.drought.chuk = 0,
  beta.wintsev.chuk = 0,
  beta.rabbit.chuk = 0,
  # alpha.sg =  rep(0,2),
  # beta.drought.sg = rep(0, 2),
  # beta.wintsev.sg = rep(0, 2),
  # beta.rabbit.sg = rep(0, 2),
  sigma.chuk = rep(1,13),
  chuk.eps = matrix(0, nrow = 13, ncol = ncol(chukar_na)-1),
  mod.chuk = rep(1,13)
)

inits <- initsFunction()

### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c("beta.spl.hunt", "mu.hunt", "pred.spl.hunt", "hunt.eps", "H", "Sigma.hunt",
                      "mu.harv",
                      'log.r.sg','C.sg',
                      "log.r.chuk", "C.chuk", "rate.chuk"))
model_test$initializeInfo()
model_test$calculate()

#Set Monitors
# pars1 <- c("alpha.hunt",
#              "Q.hunt",
#              "sig.hunt", 
#              "mu.drought.hunt",
#              "sig.drought.hunt",
#              "mu.wintsev.hunt",
#              "sig.wintsev.hunt",
#              "mu.jobs",
#              "sig.jobs",
#              "mu.income",
#              "sig.income",
#              "mu.license",
#              "sig.license",
#              "alpha.hunt",
#              "beta.drought.hunt",
#              "beta.wintsev.hunt",
#              "beta.jobs",
#              "beta.income",
#              "beta.license",
#              "beta.spl.hunt",
#              "sig.spl.hunt", 
#              "sig.une",
#              "mu.pdi",
#              "sig.gas",
#              "sig.pdi",
#              "alpha.pdi",
#              "beta.t.pdi",
#              "ar1",
#              "rho.hunt",
#              "H")
pars1 <- c(### Hunter Effort
           "alpha.hunt",
           "beta.drought.hunt",
           "beta.wintsev.hunt",
           "beta.jobs",
           "beta.income",
           "beta.license",
           "beta.spl.hunt",
           
           ### Total Harvest
           "alpha.harv",
           "beta.drought.harv",
           "beta.wintsev.harv",
           "beta.spl.harv",
           "beta.rabbit",
           
           ### Sage Grouse Wing-Bee
           "alpha.sg",
           "beta.wintsev.sg",
           "beta.drought.sg",
           "beta.rabbit.sg",
           "sg.eps",
           "mod.sg",
           
           ### Chukar Site Abundance
           "alpha.chuk",
           "beta.wintsev.chuk",
           "beta.drought.chuk",
           "beta.rabbit.chuk",
           "chuk.eps",
           "mod.chuk"
           )

# Parallel Processing Setup
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
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
  
  model_test$simulate(c("beta.spl.hunt", "mu.hunt", "pred.spl.hunt", "hunt.eps", "H", "Sigma.hunt",
                        "mu.harv",
                        'log.r.sg','C.sg',
                        "log.r.chuk", "C.chuk", "rate.chuk"))
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
save(files, file = 'model_output_FullModel.rdata')


### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList2, filename = "Traceplots - Full Model MCMC.pdf")

