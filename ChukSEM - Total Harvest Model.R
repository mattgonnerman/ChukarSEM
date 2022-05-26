### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2017 #Only need to change this to adjust the number of years

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- 2017-cutoff.y
source("./ChukarSEM - 1 Data Prep - Predict.R")

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Total Harvest ###
  mu.drought.harv ~ dnorm(0, 0.01)
  sig.drought.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.wintsev.harv ~ dnorm(0, 0.01)
  sig.wintsev.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  # mu.rabbit ~ dnorm(0, 0.01)
  # sig.rabbit ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.raven ~ dnorm(0, 0.01)
  sig.raven ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.nharrier ~ dnorm(0, 0.01)
  sig.nharrier ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
    beta.drought.harv[s] ~ dnorm(mu.drought.harv, sd = sig.drought.harv)
    beta.wintsev.harv[s] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
    # beta.rabbit.harv[s] ~ dnorm(mu.rabbit, sd  = sig.rabbit)
    beta.raven.harv[s] ~ dnorm(mu.raven, sd  = sig.raven)
    beta.nharrier.harv[s] ~ dnorm(mu.nharrier, sd  = sig.nharrier)
  } #s
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.harv[s,r] ~ dnorm(0, sd = 1)
      for(k in 1:K){
        beta.spl.harv[s,r,k] ~ dnorm(0, sd = sig.spl.harv[s,r])
      } #k
      sig.spl.harv[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Total harvest, Year 1
      
      for(t in 2:(n.year)){
        mu.harv[s,t-1,r] <- alpha.harv[s,r] +#regression formula
          beta.drought.harv[s] * pdsi[t-1,r] + #previous breeding season drought index
          beta.wintsev.harv[s] * awssi[r,t-1] + #concurrent winter severity
          # beta.rabbit.harv[s] * rabbits[t-1,r] + #concurrent rabbit harvest
          beta.raven.harv[s] * raven[t] + #prior BBS index 
          beta.nharrier.harv[s] * nharrier[t] + #prior BBS index 
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

data <- list(
  ### Total Harvest
  n.harv = upland,
  awssi = awssi, #winter severity index, scaled
  pdsi = pdsi, #Previous breeding season drought index
  # rabbits = rabbits, #Number of rabbits harvested 
  raven = as.vector(bbs.df$raven), #bbs bayes index for ravens, t = 1 is 1975
  nharrier = as.vector(bbs.df$nharrier), #bbs bayes index for northern harriers, t = 1 is 1975
  Z.harv = ZZ #Spline
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

# Wrapper Function
initsFunction <- function() list(   
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
  # N = Ni,
  beta.drought.harv = rep(0, 7),
  beta.wintsev.harv = rep(0, 7),
  # beta.rabbit.harv = rep(0, 7),
  beta.raven.harv = rep(0, 7),
  beta.nharrier.harv = rep(0, 7),
  beta.spl.harv = array(0, dim = c(7,2,12)),
  sig.spl.harv = matrix(1, ncol = 2, nrow = 7),
  mu.drought.harv = 0,
  sig.drought.harv = 1,
  mu.wintsev.harv = 0,
  sig.wintsev.harv = 1,
  # mu.rabbit = 0,
  # sig.rabbit = 1,
  mu.raven = 0,
  sig.raven = 1,
  mu.nharrier = 0,
  sig.nharrier = 1
)

inits <- initsFunction()


### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits )

model_test$simulate(c('mu.harv', 'N'))
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
  # "beta.rabbit.harv",
  "beta.raven.harv",
  "beta.nharrier.harv",
  "pred.spl.harv",
  "log.r.harv",
  "rho.harv"
  
)

### Parallel Processing Code
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
  
  model_test$simulate(c("mu.harv"))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  pars1)
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
save(files, file = 'model_output_TotHarv.rdata')

### Traceplots
# colnames(mcmcList2$chain1)
#Individual parameters
# MCMCtrace(mcmcList2, params = "alpha.hunt", plot = T, pdf = F)
#Output full pdf with all trace plots
MCMCtrace(mcmcList2, filename = "Traceplots - Total Harvest MCMC.pdf")


# harvest correlation
rho.harv.est <- MCMCsummary(mcmcList2, 'rho.harv') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList2, 'rho.harv'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 7, ncol = 7))
rho.harv.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.harv.est)){
  rho.harv.list[[rho.harv.est$Region[i]]][rho.harv.est$Species1[i], rho.harv.est$Species2[i]] <- rho.harv.est$Estimate[i]
}

rho.harv.e <- rho.harv.list[[1]]
rho.harv.w <- rho.harv.list[[2]]

rho.harv.e.plot <- ggcorrplot::ggcorrplot(rho.harv.e, lab = T) +
  labs(title = "harver Correlation - East")
ggsave(rho.harv.e.plot, filename = "CheckPlot - rho harv East.jpg", dpi = 300)
rho.harv.w.plot <- ggcorrplot::ggcorrplot(rho.harv.w, lab = T) +
  labs(title = "harver Correlation - West")
ggsave(rho.harv.w.plot, filename = "CheckPlot - rho harv West.jpg", dpi = 300)

#Comparison of Betas
require(stringr)
if(drop.rabbit == "Y"){
  check.species <- species[-c(4,7,8)]
}else{
  check.species <- species[-c(4,8)]
}
betas.vec <- c("beta.drought.harv", "beta.nharrier.harv", "beta.raven.harv", "beta.wintsev.harv")
beta.names <- c("Drought", "Harrier", "Raven", "Winter")
for(i in 1:length(betas.vec)){
  test <- MCMCsummary(mcmcList2, betas.vec[i]) %>%
    mutate(RowID = rownames(MCMCsummary(mcmcList2, betas.vec[i]))) %>%
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
  ggtitle("Effect Size - Total Harvest") +
  labs(color = "Species") + 
  theme(legend.title.align=0.5)

ggsave(hunt.betas.plot, filename = "Total Harvest - Betas - Solo.jpg", dpi = 300, width = 20, height = 20)

