### Prepare Final Data Inputs
## Splines
require(splines)
### Function that Constructs B-Spline Base
### from? = https://github.com/andrewcparnell/jags_examples/blob/master/R%20Code/jags_spline.R
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
BM1 <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z1  <- array(NA, dim = c(cut+4,nseg+2,7,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))
time <- 1:(cut+4)
for(i in 1:7){
  for(j in 1:2){
    BM1[,,i,j] <- bs_bbase(time, nseg = 10)
    Z1[,,i,j] <-  BM1[,,i,j]%*% Q1
    
  }
}

ZZ1 <- Z1
ZZ1[is.na(ZZ1)] <- 0


# Consolidate individual data objects into list for NIMBLE
data <- list(une = c(une,NA), #unemployment
             PDI = PDI, #personal disposable income
             GAS = GAS, #gas prices
             n.hunt = abind(hunters,array(NA, dim = c(7,4,2)) ,along = 2), #number of hunters (array = section)
             ZZ = ZZ1, #Spline...things...?
             res= scale(res)[,1], #residential license sailes
             wpdsi = data.matrix(abind(wpdsi,matrix(NA,1,2), along = 1))) #drought metric


### Define Constants
require(LaplacesDemon)
n.species<- dim(hunters)[1]
sig = rgamma(n.species,1,1)
Lambda = diag(sig)
I = diag(n.species) #identity matrix

# Package as list for NIMBLE
constants <- list(      n.species = 7, #number of species
                        K= 12, #number of knots for spline
                        n.region = 2, #East versus West
                        n.year = cut+4, #Number of Years
                        mu.hunt = rep(0, 7), #mean change in harvest (0 centered)
                        era = c(rep(1,19),rep(2, 27)), #Groupings for change in gas prices 
                        mean.H = apply(hunters, c(1,3), mean, na.rm = TRUE), # Mean Harvest
                        sd.H = apply(hunters, c(1,3), sd, na.rm = TRUE), # SD Harvest
                        I = abind(I,I,along = 3)) #Identity Matrix


### Set Initial Values
nu = n.species + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.species, n.species)
for (j in 1:n.species){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)

Hi <- hunters + 2500
Hi[,-1,] <- NA
Hi[7,10,] <- rpois(2,colMeans(hunters[7,-10,]))
Hi <- abind(Hi[,1:cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi <- abind(hunters,hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi[7,10,] <- Hi[7,10,]  

rho.hunt.init <- diag(7)
PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), PDI)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, GAS)

initsFunction <- function() list(  #Dan's
  GAS = GAS.inits,
  PDI = PDI.inits,
  sig.pdi = 1, 
  sig.gas = 1,
  beta.drought2 = matrix(0, 7, 2), 
  mu.drought2 = c(0,0), 
  sig.drought2 = c(1,1),
  beta.income = matrix(0, 7, 2), 
  mu.incom = c(0,0), 
  sig.incom= c(1,1),
  beta.jobs = matrix(0, 7, 2), 
  mu.jobs = c(0,0), 
  sig.jobs = c(1,1), 
  sig.wpdsi = c(1,1), 
  sig.une = 1,
  
  b0.pdi = 2.9, 
  bt.pdi = 0, 
  b0.gas = 0.9,  
  ar1 = 0, 
  bt.gas = c(0,0),
  Q = abind(Q,Q,along = 3),
  Sigma = abind(Sigma,Sigma,along = 3),
  P = abind(P,P,along = 3),
  Lambda = abind(Lambda,Lambda,along = 3),
  Delta = abind(Delta,Delta,along = 3),
  rho = abind(diag(n.species),diag(n.species),along = 3),
  
  H = Hi, 
  x = zi,
  sig.trend = matrix(1, ncol = 2, nrow = 7),
  beta.general = matrix(0, ncol = 2, nrow = 7), 
  mu.gen = c(0,0), 
  sig.gen = c(1,1),
  lbo1 = matrix(0, ncol = 2, nrow = 7),
  
  lambda1 = array(1, dim = c(7,cut+3,2) ),
  log.r1  = array(0, dim = c(7,cut+3,2) ),
  hunt.eps = array(rnorm(7*2*cut+4,0,0.1),dim = c(7,2,cut+4)),
  
  wpdsi = matrix(0,46,2),
  WPDSI = matrix(0,46,2),
  une = rep(0,46),
  UNE = rep(0,46)
)


inits <- initsFunction()


### Set Parameter Monitors
pars1 <- c('b0')
