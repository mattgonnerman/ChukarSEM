### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2017 #Only need to change this to adjust the number of years

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- 2017-cutoff.y
source("./ChukSEM - Data Prep - Predict.R")

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Predictors
  #Personal Disposable Income
  alpha.pdi ~ dnorm(0, sd = 100) #intercept
  beta.t.pdi ~ dnorm(0, sd = 100) #year
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1.pdi ~ dunif(-1,1) #Autoregressive parameter
  
  pdi.trend[1] <- alpha.pdi + beta.t.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    pdi.trend[t] <- alpha.pdi + beta.t.pdi * t
    mu.pdi[t] <- pdi.trend[t] + ar1.pdi * (PDI[t-1] - pdi.trend[t-1])
  } #t
  
  #Gas Prices
  for(i in 1:2){
    alpha.gas[i] ~ dnorm(0, sd = 10)
  }
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  # beta.gas[1] ~ dnorm(0, 0.01)
  # beta.gas[2] ~ dnorm(0, 0.01)
  
  #Relative Cost of Gas
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(alpha.gas[era.gas[t]], sd = sig.gas),0.000000000001,)
    rel.cost[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  } #t
  
  #Unemployment Rate
  sig.une ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
  } #t
  
  #Resident License sales
  sig.res ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(t in 1:n.year){
    res[t] ~ dnorm(0, sd = sig.res)
  }
  
  #Drought Index
  for(r in 1:n.region){
    sig.wpdsi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
    } #t
  } #r
  
  #Winter Severity
  beta.awssi ~ dnorm(0, sd = 100)
  alpha.awssi ~ dnorm(0, sd = 100)
  for(r in 1:n.region){
    sig.awssi[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(t in 1:n.year){
      awssi[r,t] ~ dnorm(alpha.awssi + beta.awssi*(era.awssi[t]-1), sd = sig.awssi[r])
    }
  }
  
  ################################################################################
  ### Hunter Effort ###
  mu.drought.hunt ~ dnorm(0, 0.01)
  sig.drought.hunt ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.wintsev.hunt ~ dnorm(0, 0.01)
  sig.wintsev.hunt ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.jobs ~ dnorm(0, 0.01)
  sig.jobs ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.income ~ dnorm(0, 0.01)
  sig.income ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.license ~ dnorm(0, 0.01)
  sig.license ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
    beta.drought.hunt[s] ~ dnorm(mu.drought.hunt, sd = sig.drought.hunt)
    beta.wintsev.hunt[s] ~ dnorm(mu.wintsev.hunt, sd  = sig.wintsev.hunt)
    beta.jobs[s] ~ dnorm(mu.jobs, sd = sig.jobs)
    beta.income[s] ~ dnorm(mu.income, sd  = sig.income)
    beta.license[s] ~ dnorm(mu.license, sd  = sig.license)
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
        mu.hunt[s,t,r] <- alpha.hunt[s,r] + #intercept
          beta.drought.hunt[s] * wpdsi[t,r] + #concurrent winter drought index
          beta.wintsev.hunt[s] * awssi[r,t] + #concurrent winter severity
          beta.jobs[s] * une[t] + #concurrent years unemployment
          beta.income[s] * rel.cost[t] + #PDI/Gas Price
          beta.license[s] * res[t] + #Hunting licences sold that season
          inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #spline smoothing
        
        pred.spl.hunt[s,r,t] <- inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #Derive spline smoothing for examination later
        
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
  ### Hunter Effort
  n.hunt = hunters, #Observed number of hunters for each species each year
  awssi = awssi, #winter severity index, scaled
  une = une, #BL Unemployment information for Nevada, scaled
  Z.hunt = ZZ1, #Spline 
  res = scale(res)[,1], #Residential license sales
  wpdsi = wpdsi, #winter drought index, scaled
  PDI = PDI, #personal disposable income
  GAS = GAS #gas prices
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
  
  ###Predictors
  era.gas = c(rep(1,length(1976:1994)),rep(2, length(1995:2017))), #Groupings for change in gas prices 
  era.awssi = c(rep(1,length(1976:1994)),rep(2, length(1995:2001)), rep(1, length(2002:2017))), #Groupings for change in gas prices 
  
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
  alpha.pdi = 2.9,
  beta.t.pdi = 0,
  sig.pdi = 1,
  ar1.pdi = 0,
  PDI = PDI.inits,
  #Gas Prices
  GAS = GAS.inits,
  alpha.gas = rep(.9,2),
  # beta.gas = 0,
  sig.gas = 1,
  #Unemployment
  sig.une = 1,
  une = une.init,
  #Resident License Sales
  sig.res = 1,
  res = res.init,
  #Drought Index
  sig.wpdsi = rep(1,2),
  sig.pdsi = rep(1,2),
  wpdsi = wpdsi.init,
  # pdsi = pdsi.init,
  #Winter Severity
  sig.awssi = rep(1,2),
  beta.awssi = 0,
  alpha.awssi = 0,
  awssi = awssi.init,
  
  sig.spl.hunt = matrix(1, ncol = 2, nrow = 7),
  mu.drought.hunt = 0,
  sig.drought.hunt = 1,
  mu.wintsev.hunt = 0,
  sig.wintsev.hunt = 1,
  mu.jobs = 0,
  sig.jobs = 1,
  mu.income = 0,
  sig.income = 1,
  mu.license = 0,
  sig.license = 1,
  
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
  beta.license = rep(0, 7)
)

inits <- initsFunction()

### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c('mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H', 'Sigma.hunt', 'lambda.hunt'))
model_test$initializeInfo()
model_test$calculate()

#Set Monitors
pars1 <- c("ar1.pdi",
           "alpha.hunt",
           "beta.drought.hunt",
           "beta.wintsev.hunt",
           "beta.jobs",
           "beta.income",
           "beta.license",
           "beta.spl.hunt",
           "rho.hunt",
           "H")

pars2 <- c(
  "alpha.pdi",
  "beta.t.pdi",
  "sig.pdi",
  "ar1.pdi",
  "alpha.gas",
  "sig.gas",
  # "beta.gas",
  "rel.cost",
  "sig.une",
  "une",
  "sig.res",
  "res",
  "sig.wpdsi",
  # "sig.pdsi",
  "sig.awssi",
  "beta.awssi",
  "alpha.awssi"
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
  
  model_test$simulate(c("beta.spl.hunt", "mu.hunt", "pred.spl.hunt", "hunt.eps", "H", "Sigma.hunt"))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors = pars1, monitors2 =  pars2) 
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 100000, niter = 200000, thin = 100, thin2 = 100)
  
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
save(files, file = 'model_output_HuntEff.rdata')


### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList2, filename = "Traceplots - Hunter Effort MCMC.pdf")

MCMCtrace(mcmcList1, filename = "Traceplots - Hunter Effort MCMC Betas.pdf")

# harvest correlation
rho.hunt.est <- MCMCsummary(mcmcList1, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 7, ncol = 7))
rho.hunt.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.hunt.est)){
  rho.hunt.list[[rho.hunt.est$Region[i]]][rho.hunt.est$Species1[i], rho.hunt.est$Species2[i]] <- rho.hunt.est$Estimate[i]
}

rho.hunt.e <- rho.hunt.list[[1]]
rho.hunt.w <- rho.hunt.list[[2]]

rho.hunt.e.plot <- ggcorrplot::ggcorrplot(rho.hunt.e, lab = T) +
  labs(title = "Hunter Correlation - East")
ggsave(rho.hunt.e.plot, filename = "CheckPlot - rho hunt East.jpg", dpi = 300)
rho.hunt.w.plot <- ggcorrplot::ggcorrplot(rho.hunt.w, lab = T) +
  labs(title = "Hunter Correlation - West")
ggsave(rho.hunt.w.plot, filename = "CheckPlot - rho hunt West.jpg", dpi = 300)

#Comparison of Betas
require(stringr)
if(drop.rabbit == "Y"){
  check.species <- species[-c(4,7,8)]
}else{
  check.species <- species[-c(4,8)]
}
betas.vec <- c("beta.drought.hunt", "beta.income", "beta.jobs", "beta.license", "beta.wintsev.hunt")
beta.names <- c("Drought", "Income", "Jobs", "License", "Winter")
for(i in 1:length(betas.vec)){
  test <- MCMCsummary(mcmcList1, betas.vec[i]) %>%
    mutate(RowID = rownames(MCMCsummary(mcmcList1, betas.vec[i]))) %>%
    mutate(Species = check.species[as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])"))]) %>%
    mutate(Beta = beta.names[i]) %>%
    dplyr::select(Species, Beta, Estimate = mean, LCL = '2.5%', UCL = '97.5%')
  if(i == 1){
    betas.df <- test
  }else{
    betas.df <- rbind(betas.df, test)
  }
}


hunt.betas.plot <- ggplot(data = betas.df, aes(y = Beta, x = Estimate, color = Species)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2, size = 1.5) +
  geom_point(size = 8,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = LCL, xmax = UCL),
                width = 0, size = 2,
                position = position_dodge(width = .4)) +
  theme_classic(base_size = 50) + 
  xlab("Coefficient Estimate") +
  ylab("") +
  ggtitle("Effect Size - Hunter Effort") +
  labs(color = "Species") + 
  theme(legend.title.align=0.5)

ggsave(hunt.betas.plot, filename = "HuntEff - Betas - Solo.jpg", dpi = 300, width = 20, height = 20)
