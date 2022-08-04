### Load Packages and Set Data Subset Info
lapply(c("parallel", "coda", "MCMCvis", "dplyr", "ggplot2",
         "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "LaplacesDemon"), require, character.only = T)
min.year.data <- 2002
max.year.data <- 2017
max.year.est <- 2018
cut <- years.total <- length(min.year.data:max.year.est)

###############################################
### Data Prep
## Harvest Data
chukharv_data <- read.csv('./Chukar Only/ChukarHarvestData.csv') %>%
  mutate(N = round(N),
         H = round(H)) %>%
  filter(Year %in% min.year.data:max.year.data) %>%
  arrange(County) 

county_order <- unique(chukharv_data$County)
# Regions, 1 == West, 2 == East, 3 == South
county_reg <- c(1, 1, 3, 1, 2, 3, 2, 1, 2, 3, 1, 1, 3, 1, 1, 1, 2)
n_county <- length(county_order)


chukharv_data <- rbind(chukharv_data,
                       data.frame(expand.grid(Year = (max.year.data+1):max.year.est,
                                County = county_order, N = NA, H = NA)))


chuk_hunt <- chukharv_data %>%
  arrange(County, Year) %>%
  dplyr::select(-N) %>%
  pivot_wider(values_from = H, names_from = Year) %>%
  mutate(C_Order = which(county_order == County)) %>%
  arrange(C_Order) %>%
  dplyr::select(-County, -C_Order) %>%
  as.data.frame() %>%
  mutate_all(as.numeric)

chuk_harv <- chukharv_data %>%
  arrange(County, Year) %>%
  dplyr::select(-H) %>%
  pivot_wider(values_from = N, names_from = Year) %>%
  mutate(C_Order = which(county_order == County)) %>%
  arrange(C_Order)

n.year.harv <- ncol(chuk_hunt)

## Chukar Site Abundance
chuk_site_ID <- read.csv("./Data/Chukar_Surveys_locations.csv") %>%
  dplyr::select(Survey.Location, County = COUNTY) %>%
  mutate(CountyID = match(County, county_order, nomatch = NA)) %>%
  mutate(RegionID = county_reg[CountyID])

chuksiteabun_data <- read.csv('./Data/Chukar_Surveys_data.csv') %>%
  gather("Population", "Count", 2:14)

chuksiteabun_data <- rbind(chuksiteabun_data,
                           expand.grid(Year = c(2002:2007,(max.year.data+1):max.year.est),
                                       Population = unique(chuksiteabun_data$Population),
                                       Count = NA)) %>%
  pivot_wider(values_from = Count, names_from = Year) %>%
  dplyr::rename(Survey.Location = Population) %>%
  merge(., chuk_site_ID, by = "Survey.Location", all.x = T)

survey.abun <- as.matrix(chuksiteabun_data %>% dplyr::select(as.character(1990:max.year.est)))
survey.county <- chuksiteabun_data$CountyID
survey.reg <- chuksiteabun_data$RegionID
n.year.surv <- ncol(survey.abun)

## BBS Data
bbs.df <- read.csv("./Data/bbs_indices.csv") %>%
  mutate(raven = scale(raven)[,1],
         rthawk = scale(rthawk)[,1],
         nharrier = scale(nharrier)[,1],
         pfalcon = scale(pfalcon)[,1]
  ) %>%
  filter(Year %in% min.year.data:max.year.data) %>%
  dplyr::select(Year, raven) %>%
  rbind(.,
        data.frame(expand.grid(Year = (max.year.data+1):max.year.est,
                               raven = NA)))


#Winter Severity (AWSSI)
awssi <- read.csv("./Data/Nevada AWSSI.csv") %>%
  select(Station, Region, Start, End, AWSSI) %>%
  mutate(Start = as.Date(Start, "%m/%d/%Y")) %>%
  mutate(Year = lubridate::year(Start)) %>%
  group_by(Region, Year) %>%
  filter(!is.na(AWSSI)) %>%
  summarise(AWSSI = mean(AWSSI)) %>%
  filter(Year %in% (min.year.data-1):max.year.data) %>%
  mutate(AWSSI = scale(AWSSI)[,1]) %>%
  pivot_wider(names_from = Region, values_from = AWSSI) %>%
  rbind(.,
        data.frame(expand.grid(Year = (max.year.data+1):max.year.est,
                               Eastern = NA, Southern = NA, Western = NA))) %>%
  dplyr::select(Western, Eastern, Southern) %>%
  as.matrix()


## Unemployment data
## 1976 - 2021 (Bureau of Labor, https://www.bls.gov/regions/west/nevada.htm#eag)
unemployment <- read.csv("./Data/NEvadaUnemploymentUSBL.csv", colClasses = c("integer", "character", rep("numeric", 2), rep("integer", 3), "numeric")) %>%
  filter(Period == "Sep") %>%
  filter(Year >= min.year.data) %>%
  filter(Year <= max.year.data) %>%
  mutate(Rate = unemployment/labor.force) %>%
  dplyr::select(Year, Rate) %>%
  rbind(.,
        data.frame(expand.grid(Year = (max.year.data+1):max.year.est,
                               Rate = NA)))


#################################################
### Hunter Participation Model
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
time <- 1:n.year.harv

BM1 <- array(NA, dim = c(n.year.harv,nseg+3,n_county))
Z1  <- array(NA, dim = c(n.year.harv,nseg+2,n_county))
D1 <- diff(diag(ncol(BM1[,,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))

for(i in 1:n_county){
    BM1[,,i] <- bs_bbase(time, nseg = 10)
    Z1[,,i] <-  BM1[,,i]%*% Q1
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0

data <- list(
  ### Covariates
  une = unemployment$Rate, #BL Unemployment information for Nevada, scaled
  awssi = awssi, #winter severity index, scaled
  raven = bbs.df$raven, #bbs bayes index for ravens, t = 1 is 1975
  
  ### Hunter Effort
  n.hunt = chuk_hunt, #Observed number of hunters for each species each year
  Z.hunt = ZZ1 #Spline
)

### Specify Constants
sig = rgamma(n_county,1,1)
Lambda = diag(sig)
I = diag(n_county) #identity matrix

constants <- list(
  n.region = 3,
  n.counties = n_county,
  n.year = years.total,
  n.year.harv = n.year.harv,
  K = 12,
  
  ### Predictors
  # era.awssi = c(rep(0,length(min.year.data:1994)),rep(1, length(1995:2001)), rep(0, length(2002:(max.year.est+1)))), #Groupings for change in gas prices
  
  ### Hunter Effort
  I.hunt = I,
  harv_reg_id = county_reg,
  
  ### Total harvest
  I.harv = I,
  
  ### Chukar Site Abundance
  n.site = nrow(chukar),
  n.year.chuk = ncol(chukar),
  count.surv = chuk.reg
)

### Specify Initial Values
## Hunter Effort
nu = n_county + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n_county, n_county)
for (j in 1:n_county){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)


n.hunt.i <- chuk_hunt

for(i in 1:nrow(chuk_hunt)){for(j in 1:ncol(chuk_hunt)){
  n.hunt.i[i,j] <- ifelse(is.na(chuk_hunt[i,j]), floor(mean(as.matrix(chuk_hunt[i,]), na.rm = T)), NA)
}}

### Predictors
une.init <- ifelse(is.na(unemployment$Rate), 0, NA)
awssi.init <- awssi
for(i in 1:nrow(awssi)){
  for(j in 1:3){
    awssi.init[i,j] <- ifelse(is.na(awssi[i,j]), 0, NA)
  }
}
raven.init <- as.vector(ifelse(is.na(bbs.df$raven), 0, NA))

# Wrapper Function
initsFunction <- function() list(
  ### Predictors
  #Winter Severity
  sig.awssi = rep(1,3),
  beta.awssi = rep(0,3),
  alpha.awssi = rep(0,3),
  awssi = awssi.init,
  # #Ravens
  alpha.rav = 0,
  beta.t.rav = 0,
  sig.rav = 1,
  ar1.rav = 0,
  raven = raven.init,
  #Unemployment
  sig.une = 1,
  une = une.init,
  
  
  ### Hunter Effort
  alpha.hunt = rep(1, n_county),
  beta.jobs = 0,
  beta.wintsev.hunt = 0,
  beta.spl.hunt = matrix(0, nrow = n_county, ncol = 12),
  sig.spl.hunt = rep(1, n_county),
  
  n.hunt = n.hunt.i,
  Q.hunt = Q,
  P.hunt = P,
  Lambda.hunt = diag(n_county),
  Delta.hunt = Delta,
  rho.hunt = diag(n_county)
)

inits <- initsFunction()

code <- nimbleCode( {
  ################################################################################
  ### Predictors ###
  #Ravens (Highly correlated with Prairie Falcon and RTHawk)
  alpha.rav ~ dnorm(0, 0.001) #intercept
  beta.t.rav ~ dnorm(0, 0.01) #year
  sig.rav~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1.rav ~ dunif(-1,1) #Autoregressive parameter
  
  rav.trend[1] <- alpha.rav + beta.t.rav * 1
  mu.rav[1] <- rav.trend[1]
  for(t in 2:n.year){
    rav.trend[t] <- alpha.rav + beta.t.rav * t
    mu.rav[t] <- rav.trend[t] + ar1.rav * (raven[t-1] - rav.trend[t-1])
  } #t
  
  for(t in 1:n.year){
    raven[t] ~ dnorm(mu.rav[t], sd = sig.rav)
  }
  
  #Unemployment (Jobs)
  sig.une ~ T(dt(0, pow(2.5, -2), 1), 0, )
  for (t in 1:(n.year)) {
    une[t] ~ dnorm(0, sd = sig.une)
  }
  
  #Winter Severity
  for (r in 1:3) {
    # beta.awssi[r] ~ dnorm(0, sd = 10)
    alpha.awssi[r] ~ dnorm(0, sd = 10)
    sig.awssi[r] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for (t in 1:(n.year+1)) {
      # awssi[t,r] ~ dnorm(alpha.awssi[r] + beta.awssi[r] * era.awssi[t], sd = sig.awssi[r])
      awssi[t,r] ~ dnorm(alpha.awssi[r], sd = sig.awssi[r])
    }
  }
  
  ################################################################################
  ### Hunter Effort ###
    beta.wintsev.hunt ~ dnorm(0, sd = 100)
    beta.jobs ~ dnorm(0, sd = 100)

    for(s in 1:n.counties){
      alpha.hunt[s] ~ dnorm(5, sd = 12) #sd = 1)
      sig.spl.hunt[s] ~ T(dt(0, pow(2.5,-2), 1),0,)
      for(k in 1:K){
        beta.spl.hunt[s,k] ~ dnorm(0, sd = sig.spl.hunt[s])
      } #k
      
      for(t in 1:n.year.harv){ 
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,t] <- alpha.hunt[s] + #intercept
          beta.wintsev.hunt * awssi[t + 1, harv_reg_id[s]] + 
          beta.jobs * une[t] + 
          inprod(beta.spl.hunt[s,1:K], Z.hunt[t,1:K,s]) #spline smoothing
        
        H[s,t] <- exp(hunt.eps[s,t]) #Log Link
        
        n.hunt[s,t] ~ dpois(H[s,t]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year.harv-1)){
        lambda.hunt[s,t] <- H[s,t+1]/H[s,t]
        log.r.hunt[s,t] <- log(lambda.hunt[s,t])
      }
    } #s
    
    for(t in 1:n.year.harv){
      hunt.eps[1:n.counties,t] ~ dmnorm(mu.hunt[1:n.counties,t], cov =  Sigma.hunt[1:n.counties,1:n.counties] )
    }
    
    # Correlation Matrices
    Q.hunt[1:n.counties,1:n.counties] ~ dinvwish(S = I.hunt[1:n.counties,1:n.counties], df = n.counties + 1)
    
    for(s in 1:n.counties){
      Lambda.hunt[s,s] ~ dgamma(1,1)
      Delta.hunt[s,s] <- pow(Q.hunt[s,s], -0.5)
    } #s
    
    for (s1 in 2:n.counties){
      for (s2 in 1:(s1-1)){
        Lambda.hunt[s1,s2] <- 0
        Delta.hunt[s1,s2] <- 0
      } #s2
    } #s1
    
    Sigma.hunt[1:n.counties,1:n.counties] <- Lambda.hunt[1:n.counties,1:n.counties] %*% P.hunt[1:n.counties,1:n.counties] %*% Lambda.hunt[1:n.counties,1:n.counties]  
    P.hunt[1:n.counties,1:n.counties] <- Delta.hunt[1:n.counties,1:n.counties] %*% Q.hunt[1:n.counties,1:n.counties] %*% Delta.hunt[1:n.counties,1:n.counties]
    
    for (s1 in 1:n.counties){
      for (s2 in 1:n.counties){
        rho.hunt[s1,s2] <- Sigma.hunt[s1,s2]/sqrt(Sigma.hunt[s1,s1] * Sigma.hunt[s2,s2])   
      } #s2
    } #s1
})


#Check Model Before Running
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)
model_test$simulate(c(
  'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H',
  'Sigma.hunt', 'lambda.hunt', 'log.r.hunt'
))
model_test$initializeInfo() #Want all variables initialized
model_test$calculate() #Looking for non NA value


### Run Hunter Effort Model
#Set Monitors
pars1 <- c("H")
# Parallel Processing Setup
rm(out)
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
  
  model_test$simulate(c(
    'mu.hunt', 'beta.spl.hunt', 'pred.spl.hunt', 'hunt.eps', 'H',
    'Sigma.hunt', 'lambda.hunt', 'log.r.hunt'
    ))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors = pars1) 
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 45000, niter = 50000, thin = 5, thin2 = 5)
  
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)

#Find model runtime
end_time <- Sys.time()
end_time - start_time

samples1    <- list(chain1 =  out[[1]], 
                    chain2 =  out[[2]], 
                    chain3 =  out[[3]])

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))



#Save Outputs as file
files <- list(mcmcList1,code)
save(files, file = './Chukar Only/Output - hunteffonly.rdata')



################################################
### Full Model
hunter.prime   <- MCMCpstr(mcmcList1, 'H')$H #Extract hunter numbers from Model1

nseg <- 10 #Number of spline segments
BM <- array(NA, dim = c(n.year.harv,nseg+3,17))
Z  <- array(NA, dim = c(n.year.harv,nseg+2,17))
D <- diff(diag(ncol(BM[,,1])), diff = 1)
Q <- t(D) %*% solve(D %*% t(D))

for(i in 1:17){
    BM[,,i] <- bs_bbase(hunter.prime[i,], nseg = 10)
    Z[,,i] <-  BM[,,i]%*% Q
}

ZZ <- Z
ZZ[is.na(ZZ)] <- 0


data <- list(
  ### Covariates
  une = unemployment$Rate, #BL Unemployment information for Nevada, scaled
  awssi = awssi, #winter severity index, scaled
  raven = bbs.df$raven, #bbs bayes index for ravens, t = 1 is 1975
  
  ### Hunter Effort
  n.hunt = chuk_hunt, #Observed number of hunters for each species each year
  Z.hunt = ZZ1, #Spline
  
  ### Total Harvest
  n.harv = chuk_harv,
  Z.harv = ZZ, #Spline
  
  ### Chukar Site Abundance
  n.chuk = survey.abun
)


### Specify Constants
n.county<- dim(chuk_hunt)[1]
sig = rgamma(n.county,1,1)
Lambda <- Lambda2 <- diag(sig)
I2 <- I <- diag(n.county) #identity matrix

constants <- list(
  n.region = 2,
  n.species = nrow(chuk_harv),
  n.year = ncol(chuk_hunt),
  K = 12,
  
  ### Predictors
  era.awssi = c(rep(0,length(1975:1994)),rep(1, length(1995:2001)), rep(0, length(2002:final.y))), #Groupings for change in gas prices
  
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

econ.inits <- econ_data %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))
bbs.inits <- as.data.frame(bbs.df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))
awssi.inits <- as.matrix(as.data.frame(awssi.df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))
pdsi.inits <- as.matrix(as.data.frame(pdsi_df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))

# Wrapper Function
initsFunction <- function() list(
  ### Covariates
  gas = econ.inits$Gas.May,
  une = econ.inits$Une,
  res = econ.inits$Licenses,
  pdi = econ.inits$PDI,
  ravens = bbs.inits$raven,
  rthawk = bbs.inits$rthawk,
  nharr = bbs.inits$nharrier,
  pfal = bbs.inits$pfalcon,
  awssi = awssi.inits,
  pdsi = pdsi.inits,
  sig.econ.pred = rep(1,4),
  mu.econ.p = rep(0,4),
  zeta.econ = rep(1,ncol(econ_data)),
  pred.econ.prime = rep(0, cut),
  mu.econ.pred = matrix(0, 4, cut),
  sig.bbs.pred = rep(1,4),
  mu.bbs.p = rep(0,4),
  zeta.bbs = rep(1,ncol(bbs.df)),
  pred.bbs.prime = rep(0, cut),
  mu.bbs.pred = matrix(0, 4, cut),
  beta.awssi = 0,
  alpha.awssi = 0,
  sig.awssi = rep(1,2),
  sig.drought = rep(1,2),
  
  ### Covariates
  mu.econ = 0,
  sig.econ = 1,
  beta.econ.hunt = rep(0, 5),
  beta.spl.hunt = array(0, dim = c(5,2,12)),
  sig.spl.hunt = matrix(1, ncol = 2, nrow = 5),
  
  mu.wintsev.harv = 0,
  sig.wintsev.harv = 1,
  beta.wintsev.harv = rep(0, 5),
  mu.bbs = 0,
  sig.bbs = 1,
  beta.bbs.harv = rep(0, 5),
  mu.pdsi = 0,
  sig.pdsi = 1,
  beta.pdsi.harv = rep(0, 5),
  mu.hunter.harv = 0,
  sig.hunter.harv = 1,
  beta.hunter.harv = matrix(0, 5, 2),
  # beta.spl.harv = array(0, dim = c(5,2,12)),
  # sig.spl.harv = matrix(1, ncol = 2, nrow = 5),
  
  ### Hunter Effort
  n.hunt = n.hunt.i,
  Q.hunt = abind(Q,Q,along = 3),
  P.hunt = abind(P,P,along = 3),
  Lambda.hunt = abind(diag(n.species),diag(n.species),along = 3),
  Delta.hunt = abind(Delta,Delta,along = 3),
  rho.hunt = abind(diag(n.species),diag(n.species),along = 3),
  alpha.hunt = matrix(0, ncol = 2, nrow = 5),
  sig.hunt = matrix(1, ncol = 2, nrow = 5),
  
  ### Total Harvest
  n.harv = n.harv.i,
  Q.harv = abind(Q2,Q2,along = 3),
  P.harv = abind(P2,P2,along = 3),
  Lambda.harv = abind(diag(n.species),diag(n.species),along = 3),
  Delta.harv = abind(Delta2,Delta2,along = 3),
  rho.harv = abind(diag(n.species),diag(n.species),along = 3),
  alpha.harv = matrix(0, ncol = 2, nrow = 5),
  sig.harv = matrix(1, ncol = 2, nrow = 5),
  log.r.harv = array(0, dim = c(5,(cut)-1,2) ),
  N = Ni,
  
  ### Chukar Site Abundance
  theta.chuk = rep(1,2),
  mod.chuk = rep(1,2),
  n.chuk = as.matrix(chukar_na),
  log.r.chuk = r.chuk.init
)

inits <- initsFunction()

### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Predictors ###
  # Economic Predictors
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
  
  # Predator Predictor
  for(i in 1:4){
    sig.bbs.pred[i] ~ dgamma(1,1)
    mu.bbs.p[i] ~ dnorm(0, 1)
  }
  zeta.bbs[1] <- 1
  for(i in 2:4){
    zeta.bbs[i] ~ dnorm(0, 1)
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
  
  # Winter Severity
  for (r in 1:3) {
    beta.awssi[r] ~ dnorm(0, sd = 10)
    alpha.awssi[r] ~ dnorm(0, sd = 10)
    sig.awssi[r] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for (t in 1:(n.year+1)) {
      awssi[t,r] ~ dnorm(alpha.awssi[r] + beta.awssi[r] * era.awssi[t], sd = sig.awssi[r])
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
  sig.econ ~ T(dt(0, pow(2.5, -2), 1), 0, )
  beta.econ.hunt ~ dnorm(0, sd = sig.econ)
  
  for(s in 1:n.counties){
    alpha.hunt[s] ~ dnorm(5, sd = 12) #sd = 1)
    sig.spl.hunt[s] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(k in 1:K){
      beta.spl.hunt[s,k] ~ dnorm(0, sd = sig.spl.hunt[s])
    } #k
    
    for(t in 1:n.year.harv){ 
      #Unlinked estimate of Hunter Numbers
      mu.hunt[s,t] <- alpha.hunt[s] + #intercept
        beta.econ.hunt * pred.econ.prime[t] + 
        inprod(beta.spl.hunt[s,1:K], Z.hunt[t,1:K,s]) #spline smoothing
      
      H[s,t] <- exp(hunt.eps[s,t]) #Log Link
      
      n.hunt[s,t] ~ dpois(H[s,t]) #Number of hunters follows Poisson
    } #t
    
    for(t in 1:(n.year.harv-1)){
      lambda.hunt[s,t] <- H[s,t+1]/H[s,t]
      log.r.hunt[s,t] <- log(lambda.hunt[s,t])
    }
  } #s
  
  for(t in 1:n.year.harv){
    hunt.eps[1:n.counties,t] ~ dmnorm(mu.hunt[1:n.counties,t], cov =  Sigma.hunt[1:n.counties,1:n.counties] )
  }
  
  # Correlation Matrices
  Q.hunt[1:n.counties,1:n.counties] ~ dinvwish(S = I.hunt[1:n.counties,1:n.counties], df = n.counties + 1)
  
  for(s in 1:n.counties){
    Lambda.hunt[s,s] ~ dgamma(1,1)
    Delta.hunt[s,s] <- pow(Q.hunt[s,s], -0.5)
  } #s
  
  for (s1 in 2:n.counties){
    for (s2 in 1:(s1-1)){
      Lambda.hunt[s1,s2] <- 0
      Delta.hunt[s1,s2] <- 0
    } #s2
  } #s1
  
  Sigma.hunt[1:n.counties,1:n.counties] <- Lambda.hunt[1:n.counties,1:n.counties] %*% P.hunt[1:n.counties,1:n.counties] %*% Lambda.hunt[1:n.counties,1:n.counties]  
  P.hunt[1:n.counties,1:n.counties] <- Delta.hunt[1:n.counties,1:n.counties] %*% Q.hunt[1:n.counties,1:n.counties] %*% Delta.hunt[1:n.counties,1:n.counties]
  
  for (s1 in 1:n.counties){
    for (s2 in 1:n.counties){
      rho.hunt[s1,s2] <- Sigma.hunt[s1,s2]/sqrt(Sigma.hunt[s1,s1] * Sigma.hunt[s2,s2])   
    } #s2
  } #s1
   
  
  
  ################################################################################
  ### Total Harvest ###
  beta.wintsev.harv ~ dnorm(0, sd  = 10)
  beta.bbs.harv ~ dnorm(0, sd  = 10)
  beta.pdsi.harv ~ dnorm(0, sd  = 10)

    for(s in 1:n.counties){
      alpha.harv[s] ~ dnorm(0, sd = 12)
      for(k in 1:K){
        beta.spl.harv[s,k] ~ dnorm(0, sd = sig.spl.harv[s])
      } #k
      sig.spl.harv[s] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      for(t in 1:n.year.harv){
        #Unlinked estimate of Hunter Numbers
        mu.harv[s,t] <- alpha.harv[s] +#regression formula
          beta.wintsev.harv[s] * awssi[t, harv_reg_id[s]] + #Previous winter severity (Affecting Survival)
          beta.pdsi.harv[s] * pdsi[t,harv_reg_id[s]] + # Same year, spring/summer drought (Affecting Survival/Reproduction)
          beta.bbs.harv[s] * pred.bbs.prime[t] + # Latent predator index (Affecting Reproduction)
          inprod(beta.spl.harv[s,1:K], Z.harv[t,1:K,s]) #spline smoothing
        
        N[s,t] <- exp(harv.eps[s,t]) #Log Link
        
        n.harv[s,t] ~ dpois(N[s,t]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.harv[s,t] <- N[s,t+1]/N[s,t]
        log.r.harv[s,t] <- log(lambda.harv[s,t])
      }
    } #s
    
    for(t in 1:n.year.harv){
      harv.eps[1:n.counties,t] ~ dmnorm(mu.harv[1:n.counties,t], cov =  Sigma.harv[1:n.counties,1:n.counties] )
    }
    
    ### Correlation Matrices
    Q.harv[1:n.counties,1:n.counties] ~ dinvwish(S = I.harv[1:n.counties,1:n.counties], df = n.counties + 1)
    
    for(s in 1:n.counties){
      Lambda.harv[s,s] ~ dgamma(1,1)
      Delta.harv[s,s] <- pow(Q.harv[s,s], -0.5)
    } #s
    
    for (s1 in 2:n.counties){
      for (s2 in 1:(s1-1)){
        Lambda.harv[s1,s2] <- 0
        Delta.harv[s1,s2] <- 0
      } #s2
    } #s1
    
    Sigma.harv[1:n.counties,1:n.counties] <- Lambda.harv[1:n.counties,1:n.counties] %*% P.harv[1:n.counties,1:n.counties] %*% Lambda.harv[1:n.counties,1:n.counties]  
    P.harv[1:n.counties,1:n.counties] <- Delta.harv[1:n.counties,1:n.counties] %*% Q.harv[1:n.counties,1:n.counties] %*% Delta.harv[1:n.counties,1:n.counties]  
    
    for (s1 in 1:n.counties){
      for (s2 in 1:n.counties){
        rho.harv[s1,s2] <- Sigma.harv[s1,s2,r]/sqrt(Sigma.harv[s1,s1] * Sigma.harv[s2,s2])   
      } #s2
    } #s1
  
  ################################################################################
  ### Chukar Site Abundance ###
  for(p in 1:n.site){
    theta.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.surv){
      log.r.chuk[p,t-1] <-  log.r.harv[count.surv[p], t] #log.r.harv[t=13] is 1990-1991
      
      C.chuk[p,t] <- exp(log.r.chuk[p,t-1]) * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[p]/(theta.chuk[p] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(prob = rate.chuk[p,t-1], size = theta.chuk[p]) #obs. # of chukars follow neg-bin
    } #t
  } #p
  
  ################################################################################
  ### Birds per Hunter
  for(t in 1:n.year.harv){
    for(s in 1:n.counties){
        BPH[s,t] <- n.harv[s,t]/n.hunt[s,t]
      }
    }
})