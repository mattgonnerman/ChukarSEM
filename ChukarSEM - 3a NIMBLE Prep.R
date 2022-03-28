### Prepare Final Data Inputs
## Splines
# Hunter Numbers
hunter.prime   <- MCMCpstr(mcmcList2, 'H')$H #Extract hunter numbers from Model1

nseg <- 10 #Number of spline segments
BM <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z  <- array(NA, dim = c(cut+4,nseg+2,7,2))
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
time <- 1:(cut+4)

BM1 <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z1  <- array(NA, dim = c(cut+4,nseg+2,7,2))
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


###
chukar <- df_ch[,-c(1:19)]

library(LaplacesDemon)

n.species<- dim(hunters)[1]

sig = rep(1, n.species)
Lambda = diag(sig)
I = diag(n.species)
nu = n.species + 1
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.species, n.species)

for (j in 1:n.species){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma = Lambda %*% P %*% Lambda

sig2 = rgamma(n.species, 1, 1)
Lambda2 = diag(sig2)
I2 = diag(n.species)
nu = n.species + 1
Q2 = rinvwishart(nu, I)    # note output is an array
Delta2 = matrix(0, n.species, n.species)

for (j in 1:n.species){
  Delta2[j,j] = Q2[j,j]^(-0.5)
}
P2 = Delta2 %*% Q2 %*% Delta2
Sigma2 = Lambda2 %*% P2 %*% Lambda2

### Package Data for model
data <- list(n.harv = abind(upland,array(NA, dim = c(7,4,2)) ,along = 2), #Number harvested
             une = c(une,NA),
             n.hunter = abind(hunters,array(NA, dim = c(7,4,2)) ,along = 2), 
             Z = ZZ, 
             ZZ = ZZ1, 
             pdsi = data.matrix(abind(pdsi,matrix(NA,2,2), along = 1)),
             wpdsi = data.matrix(abind(wpdsi,matrix(NA,1,2), along = 1)), 
             n.chuk =  data.matrix(chukar))


### Specify Constants
constants <- list(n.species = 7,
                  K= 12,
                  n.region = 2,
                  n.year = cut+4,
                  mu.hunt = rep(0, 7),
                  mu.bird = rep(0, 7),
                  n.site = 13,
                  diff = 14,
                  n.yr = ncol(chukar),
                  mean.H = apply(hunter.prime, c(1,3), mean, na.rm = TRUE),
                  sd.H = apply(hunter.prime, c(1,3), sd, na.rm = TRUE),
                  I = abind(I,I,along = 3), 
                  I2 =abind(I2,I2,along = 3)
                  )


### Specify Initial values
Hi <- hunters + 2500
Hi[,-1,] <- NA
Hi[7,10,] <- rpois(2,colMeans(hunters[7,-10,]))
Hi <- abind(Hi[,1:cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi <- abind(hunters,hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi[7,10,] <- Hi[7,10,]  

Ni <- upland + 5000
Ni[,-1,] <- NA
Ni[7,10,] <- rpois(2,colMeans(upland[7,-10,]))
Ni <- abind(Ni[,1:cut,],upland[,cut,],upland[,cut,],upland[,cut,],upland[,cut,], along = 2)
yi <- abind(upland,upland[,cut,],upland[,cut,],upland[,cut,],upland[,cut,], along = 2)
yi[7,10,] <- Ni[7,10,]  

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
bird.inits <- array(0, dim = c(7,2,43))
bird.inits[,,1:14] <- NA

rho.hunt.init <- diag(7)

# Package Initial Values for NIMBLE
initsFunction <- function() list( 
  phi = matrix(rbeta(14,1,1),7,2),
  beta.drought = matrix(0, 7, 2), 
  mu.drought = c(0,0), 
  sig.drought = c(1,1),
  beta.drought2 = matrix(0, 7, 2), 
  mu.drought2 = c(0,0), 
  sig.drought2 = c(1,1),
  beta.jobs = matrix(0, 7, 2), 
  mu.jobs = c(0,0),
  sig.jobs = c(1,1),
  
  C = chukar_na + 50,
  n.chuk =  chukar_na,
  sigma.chuk = rep(1,13),
  
  X0 = matrix(5, ncol = 2, nrow = 7), 
  theta2 = rep(1,13),
  
  Q = abind(Q,Q,along = 3),
  # Sigma = abind(Sigma,Sigma,along = 3),
  P = abind(P,P,along = 3),
  Lambda = abind(diag(n.species),diag(n.species),along = 3),
  Delta = abind(Delta,Delta,along = 3),
  rho = abind(diag(n.species),diag(n.species),along = 3),
  
  Q2 = abind(Q2,Q2,along = 3),
  # Sigma2 = abind(Sigma2,Sigma2,along = 3),
  P2 = abind(P2,P2,along = 3),
  Lambda2 = abind(diag(n.species),diag(n.species),along = 3),
  Delta2 = abind(Delta2,Delta2,along = 3),
  rho2 = abind(diag(n.species),diag(n.species),along = 3),
  
  bird.eps = bird.inits,
  sig.bird =  matrix(.25, ncol = 7, nrow = 2),
  N = Ni,  
  H = Hi,
  z= zi, 
  y = yi,
  sig.pressure = matrix(1, ncol = 2, nrow = 7),  sig.trend = matrix(1, ncol = 2, nrow = 7),
  lbo1 = matrix(0, ncol = 2, nrow = 7),
  lbo = matrix(0, ncol = 2, nrow = 7),
  lbo2 =  rep(0,13),
  
  lambda1 = array(1, dim = c(7,cut+3,2) ),
  log.r1  = array(0, dim = c(7,cut+3,2) ),
  hunt.eps = array(rnorm(7*2*cut+4,0,0.1),dim = c(7,2,cut+4)),
  pdsi = matrix(0,46,2), 
  wpdsi = matrix(0,46,2), 
  une = rep(0,46), 
  
  lambda = array(1, dim = c(7,cut+3,2) ),
  log.r  = array(0, dim = c(7,cut+3,2) )
)
