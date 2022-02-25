### Final Data Prep
nseg <- 10

BM1 <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z1  <- array(NA, dim = c(cut+4,nseg+2,7,2))
D1 <- diff(diag(ncol(BM1[,,1,1])), diff = 1)
Q1 <- t(D1) %*% solve(D1 %*% t(D1))

require(splines)
### UNSURE WHAT FUNCTION DOES
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
I = diag(n.species)
nu = n.species + 1
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.species, n.species)

for (j in 1:n.species){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma = Lambda %*% P %*% Lambda

Sigma = Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)

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
Hi <- hunters + 2500
Hi[,-1,] <- NA
Hi[7,10,] <- rpois(2,colMeans(hunters[7,-10,]))
Hi <- abind(Hi[,1:cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi <- abind(hunters,hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi[7,10,] <- Hi[7,10,]  

rho.hunt.init <- diag(7)
PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), PDI)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, GAS)

initsFunction <- function() list( 
  beta.drought2 = matrix(0, 7, 2), mu.drought2 = c(0,0), sig.drought2 = c(1,1),
  beta.jobs = matrix(0, 7, 2), mu.jobs = c(0,0), sig.jobs = c(1,1),
  
  X0 = matrix(5, ncol = 2, nrow = 7), 
  theta2 = rep(1,13),
  
  Q = abind(Q,Q,along = 3),
  # Sigma = abind(Sigma,Sigma,along = 3),
  P = abind(P,P,along = 3),
  Lambda = abind(diag(n.species),diag(n.species),along = 3),
  Delta = abind(Delta,Delta,along = 3),
  rho = abind(diag(n.species),diag(n.species),along = 3)
)

inits <- initsFunction()


### Set Parameter Monitors
pars1 <- c('b0')

###################################

### Unsure what function does, doesn't appear to be used anywhwere in original code
# Posdef <- function (n, ev = runif(n, 0, 10)) 
# {
#   Z <- matrix(ncol=n, rnorm(n^2))
#   decomp <- qr(Z)
#   Q <- qr.Q(decomp) 
#   R <- qr.R(decomp)
#   d <- diag(R)
#   ph <- d / abs(d)
#   O <- Q %*% diag(ph)
#   Z <- t(O) %*% diag(ev) %*% O
#   return(Z)
# }