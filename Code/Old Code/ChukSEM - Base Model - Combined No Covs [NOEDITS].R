### Run Initial Data Management
source("./ChukarSEM - 1 Data Prep.R")

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Hunter Effort ###
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.hunt[s,r] ~ dnorm(5, sd = 1)
      
      for(t in 1:(n.year)){ 
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,r,t] <- alpha.hunt[s,r] #No predictors
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
  for(r in 1:n.region){
    for(s in 1:n.species){
      # Regression Coefficients
      alpha.harv[s,r] ~ dnorm(0, sd = 1)
      
      # Process Model
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Total harvest, Year 1
      
      for(t in 2:(n.year)){
        mu.harv[s,t-1,r] <- alpha.harv[s,r] #regression formula
        
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
  ### Chukar Site Abundance ###
  for(p in 1:n.site){
    alpha.chuk[p] ~ dlogis(0,1) #Random Site Intercept
    
    theta.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod.chuk[p] ~ dlogis(0,1)
    
    for(t in 2:n.year.chuk){
      mu.chuk[p, t-1] <- mod.chuk[p] * log.r.harv[3, t+12, reg.chuk[p]] #log.r.harv[t=14] is 1990-1991
      chuk.eps[p,t-1]  ~ dnorm(mu.chuk[p, t-1], sd = sigma.chuk[p]) #Coeffient for annual change
      
      log.r.chuk[p,t-1] <- alpha.chuk[p] + chuk.eps[p,t-1] #Unlinked change in abundance, log.harv.chuk[t=1] is 1990-1991
      lambda.chuk[p,t-1] <- exp(log.r.chuk[p,t-1]) #Change in abundance
      
      C.chuk[p,t] <- lambda.chuk[p,t-1] * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[p]/(theta.chuk[p] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(rate.chuk[p,t-1], theta.chuk[p]) #obs. # of chukars follow neg-bin
    } #t
  } #p  
  
  
  ################################################################################
  ### Sage Grouse Wing-Bee ###
  ### Chukar Site Abundance
  for(r in 1:n.region){
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C.sg[r,1] ~ dpois(wing.b[r,1]) #Equivalent of Poisson lambda

    alpha.sg[r] ~ dlogis(0,1) #Intercept for
    sigma.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod.sg[r] ~ dlogis(0,1)

    for(t in 2:n.years.sg){
      mu.sg[r, t-1] <- mod.sg[r] * log.r.harv[7, t+26, r] #log.r.harv[t=28] is 2004-2005
      sg.eps[r,t-1] ~ dnorm(mu.sg[r, t-1], sd = sigma.sg[r])

      log.r.sg[r,t-1] <- alpha.sg[r] + sg.eps[r,t-1] #Variation in growth, log.harv.sg[t=1] is 2004-2005
      lambda.sg[r,t-1] <- exp(log.r.sg[r,t-1]) # Tranform between finite and instantaneous growth rate

      C.sg[r,t] <- lambda.sg[r,t-1] * C.sg[r,t-1] #Change in Count over time

      rate.sg[r,t-1] <- theta.sg[r]/(theta.sg[r] + C.sg[r,t]) #NB rate
      wing.b[r,t] ~ dnegbin(rate.sg[r,t-1], theta.sg[r]) #Wing-Bee data follows negative binomial
    } #t
  } #r
  
})

### Specify Data Inputs
data <- list(
  ### Hunter Effort
  n.hunt = hunters,
  
  ### Total Harvest
  n.harv = upland,
  
  ### Chukar Site Abundance
  n.chuk = data.matrix(chukar),
  
  ### Sage Grouse Wing-Bee
  wing.b = wing.b
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
  
  ### Hunter Effort
  I.hunt = abind(I,I,along = 3),
  
  ### Total Harvest
  I.harv = abind(I2,I2,along = 3),
  
  ### Chukar Site Abundance
  n.site = 13,
  n.year.chuk = ncol(chukar),
  reg.chuk = chuk.reg,
  
  ### Sage Grouse Wing-Bee
  n.years.sg = n.years.sg
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
  
  
  ### Total Harvest
  n.harv = n.harv.i,
  Q.harv = abind(Q2,Q2,along = 3),
  # Sigma.harv = abind(Sigma2,Sigma2,along = 3),
  P.harv = abind(P2,P2,along = 3),
  Lambda.harv = abind(diag(n.species),diag(n.species),along = 3),
  Delta.harv = abind(Delta2,Delta2,along = 3),
  rho.harv = abind(diag(n.species),diag(n.species),along = 3),
  alpha.harv = matrix(0, ncol = 2, nrow = 7),
  sig.harv = matrix(1, ncol = 2, nrow = 7),
  hunt.eps = array(rnorm(7*2*cut+4,0,0.1),dim = c(7,2,cut+4)),
  log.r.harv = array(0, dim = c(7,cut+3,2) ),
  N = Ni,
  
  ### Chukar Site Abundance
  alpha.chuk = rep(0, 13),
  theta.chuk = rep(1,13),
  sigma.chuk = rep(1,13),
  chuk.eps = matrix(0, nrow = 13, ncol = ncol(chukar_na)-1),
  C.chuk = chukar_na + 50,
  n.chuk =  n.chuk.i,
  mod.chuk = rep(1,13),
  
  ### Sage Grouse Wing-Bee
  alpha.sg =  rep(1,2),
  theta.sg = rep(1,2),
  sigma.sg = rep(1,2),
  sg.eps = matrix(0, nrow = 2, ncol = n.years.sg-1),
  C.sg = C.sg.i,
  mod.sg = rep(1,2)
  )

inits <- initsFunction()


### Check Model Code
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c("alpha.sg", "sg.eps", "C.sg", "theta.sg", "mod.sg", "rate.sg"))
model_test$initializeInfo()
model_test$calculate()


### Parallel Processing Code
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants")) #identify what is to be exported to each cluster

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
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  c("alpha.hunt",
                                                           "Q.hunt",
                                                           "sig.hunt",
                                                           
                                                           "alpha.harv",
                                                           "Q.harv",
                                                           "sig.harv",
                                                           
                                                           "alpha.chuk",
                                                           "theta.chuk",
                                                           "sigma.chuk",
                                                           "mod.chuk",
                                                           
                                                           "alpha.sg",
                                                           "theta.sg",
                                                           "sigma.sg",
                                                           "mod.sg"
                                                           )) 
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

# Traceplots
colnames(mcmcList2$chain1)
MCMCtrace(mcmcList2, params = "Q.harv", plot = T, pdf = F)
