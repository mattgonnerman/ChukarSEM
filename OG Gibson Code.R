set.seed(2)

lapply(c("dplyr", "ggplot2", "reshape2", "reshape"), require, character.only = T) 

setwd('E:/Maine Drive/Projects/Gibson Chukkar/Code/PDSI')

temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)

CC    <- myfiles[[1]]
ELKO  <- myfiles[[2]]
ELY   <- myfiles[[3]]
EUR   <- myfiles[[4]]
FAL   <- myfiles[[5]]
RL    <- myfiles[[6]]
SS    <- myfiles[[7]]
TO    <- myfiles[[8]]
WIN   <- myfiles[[9]]

CC  $Station <- rep('CC')
ELKO$Station <- rep('ELKO')
ELY $Station <- rep('ELY')
EUR $Station <- rep('EUR')
FAL $Station <- rep('FAL')
RL  $Station <- rep('RL')
SS  $Station <- rep('SS')
TO  $Station <- rep('TO')
WIN $Station <- rep('WIN')



df <- rbind(CC,ELKO,ELY,EUR,FAL,RL,SS ,TO ,WIN )

df$TimeStep <- NULL

df$Year <- format(as.Date(df$Date, format="%m/%d/%Y"),"%Y")
df$Month <- format(as.Date(df$Date, format="%m/%d/%Y"),"%m")


df$Region <- recode(df$Station, 'CC' = 'West','ELKO' = 'East','ELY' = 'East','EUR' = 'East','FAL' = 'West','RL' = 'East','SS' = 'East','TO' = 'West','WIN' = 'West')
df$Season <- recode(df$Month,`01` = 'Winter',`02` = 'Winter',`03` = 'Winter',`04` = 'Early',`05` = 'Early',`06` = 'Early',`07` = 'Late',`08` = 'Late',`09` = 'Late',`10` = 'Winter',`11` = 'Winter',`12` = 'Winter' )

df$F.Year <- ifelse(df$Month == '10' | df$Month == '11' | df$Month == '12', as.numeric(as.character(df$Year)) + 1,as.numeric(as.character(df$Year)))


df_agg <- df %>% 
  group_by(Station,F.Year,Season) %>%
  summarize(m.PDSI = mean(PDSI, na.rm = TRUE))


df_cast <- cast(df_agg, F.Year + Season ~ Station, value = 'm.PDSI')

ggplot(df_agg, aes(x = F.Year, y = m.PDSI, color = Season)) +
  geom_point()+
  facet_wrap(~Station, ncol = 3) + theme_bw()

setwd('E:/Maine Drive/Projects/Gibson Chukkar/Code')

survey_data <- read.csv('Chukar_Surveys_data.csv')
harvest_data <- read.csv('all_species_harvest_data.csv')
harvest_data$Year <- as.numeric(as.character(harvest_data$Year))
harvest_data$Animals <- as.numeric(as.character(harvest_data$Animals))
harvest_data$Hunters <- as.numeric(as.character(harvest_data$Hunters))
harvest_data$H.Days <- as.numeric(as.character(harvest_data$H.Days  ))

library(dplyr)
library(reshape2)
library(tidyverse)
library(jagsUI)
library(reshape)
library(reshape)
library(nimble)

df <- gather(survey_data, "Population", "Count", 2:14)

site_list <- unique(df$Population)
year_list <- rep(1972:2017)

df_full <- data.frame("Population" = rep(site_list, times = length(year_list)), "Year" = rep(year_list, each = length(site_list)))

df_full <- left_join(df_full, df)

df_ch <- dcast(df_full, Population ~ Year)

ggplot() +
  geom_boxplot(data = df_full, aes(x = Year, y = Count, group = Year), fill = 'dodgerblue') +
  #geom_smooth(data = df_full, aes(x = Year, y = Count), fill = 'black',alpha = .5, color = 'dodgerblue1') +
  labs(y = 'Chukar covey survey counts') +
  theme_bw()

harvest_data$Species_Code <- recode(harvest_data$Species, 
                                    `SAGE_GROUSE` = "SAGR",
                                    `SOOTY_GROUSE` = "BLGR",
                                    `RUFFED_GROUSE` = "RUGR",
                                    `CHUKAR` = "CHUK",
                                    `HUNGARIAN_PART` = 'HUPA',
                                    `CALIFORNIA_QUAIL` = 'CAQU',
                                    `MOUNTAIN_QUAIL` = 'MOQU',
                                    `PHEASANTS` = 'PHEA',
                                    `RABBIT` = 'RABB',
                                    `DOVE` = 'DOVE',
                                    `DUSKY_GROUSE` = 'BLGR',
                                    `GAMBEL'S_QUAIL` = 'GAQU'
)

harvest_datae <- subset(harvest_data, grepl('Eastern', Section))

animal_wide <- cast(harvest_datae, Species_Code  ~ Year,  value = 'Animals', fun.aggregate = 'sum', na.rm = TRUE)
animal_wide[animal_wide==0] <- NA
species <- animal_wide$Species_Code 
animal_wide$Species_Code  <- NULL
animal_wide <- data.matrix(animal_wide)

hunters_wide <- cast(harvest_datae, Species_Code  ~ Year, value = 'Hunters', fun.aggregate = 'sum', na.rm = TRUE)
hunters_wide[hunters_wide==0] <-NA
hunters_wide$Species_Code  <- NULL
hunters_wide<- data.matrix(hunters_wide)

harvest_dataw <- subset(harvest_data, grepl('Western', Section))

animal_widew <- cast(harvest_dataw, Species_Code  ~ Year,  value = 'Animals', fun.aggregate = 'sum', na.rm = TRUE)
animal_widew[animal_widew==0] <- NA
speciesw <- animal_widew$Species_Code 
animal_widew$Species_Code  <- NULL
animal_widew <- data.matrix(animal_widew)

hunters_widew <- cast(harvest_dataw, Species_Code  ~ Year, value = 'Hunters', fun.aggregate = 'sum', na.rm = TRUE)
hunters_widew[hunters_widew==0] <-NA
hunters_widew$Species_Code  <- NULL
hunters_widew<- data.matrix(hunters_widew)

animal_widew <- animal_widew[-6,]
hunters_widew <- hunters_widew[-6, ]

library(abind)
upland <- abind(animal_widew,animal_wide, along = 3)
hunters <- abind(hunters_widew, hunters_wide, along = 3)

hunters <-hunters[-c(4,8),,]
upland <- upland[-c(4,8),,]

df <- cbind.data.frame(year = c(1976:2017), chukar = c(upland[2,,1]), hunter = c(hunters[2,,1])) 

ggplot(data = df) +
  geom_point(aes(x = year, y = log(chukar)), size = 3, shape = 21, fill = 'dodgerblue') +
  geom_point(aes(x = year, y = log(hunter)), size = 3, shape = 21, fill = 'goldenrod') +
  labs(y = 'Numbers of hunters and chukar (log-scaled)', x = 'Year') +
  theme_bw()

upland[7,10,] <- NA # Season closed
hunters[7,10,] <- NA # Season closed 

eastern_pdsi <- read.csv('Eastern_PDSIZ.csv')
eastern_pdsi <- subset(eastern_pdsi, Year > 1974)
winter_pdsie <- rep(NA, times = dim(eastern_pdsi)[1] - 1)

for(i in 1:length(winter_pdsie)){
  winter_pdsie[i] <- sum(eastern_pdsi[i, 12:13], eastern_pdsi[i+1,2:4] )/5
}

ebreeding_means <- data.frame(year = eastern_pdsi$Year, breeding_drought = rowMeans(eastern_pdsi[,5:9]))

eastern_pdsi <- subset(eastern_pdsi, Year > 1975)
ebreeding_means <- data.frame(year = eastern_pdsi$Year, breeding_drought = rowMeans(eastern_pdsi[,5:9]))

western_pdsi <- read.csv('Western_PDSIZ.csv')

western_pdsi <- subset(western_pdsi, Year > 1974)
winter_pdsiw <- rep(NA, times = dim(western_pdsi)[1] - 1)

for(i in 1:length(winter_pdsiw)){
  winter_pdsiw[i] <- sum(western_pdsi[i, 12:13], western_pdsi[i+1,2:4] )/5
}

western_pdsi <- subset(western_pdsi, Year > 1975)
wbreeding_means <- data.frame(year = western_pdsi$Year, breeding_drought = rowMeans(western_pdsi[,5:9]))

cut <- 42

pdsi  <- data.frame(east = ebreeding_means$breeding_drought[1:44],west = wbreeding_means$breeding_drought[1:44] )
wpdsi <- data.frame(east = winter_pdsie[1:(44+1)], west = winter_pdsiw[1:(44+1)])

#1976 - 2021
unemployment <- c(8.825,6.766666667,4.416666667,4.866666667,6.316666667,7.316666667,10.05,9.841666667,7.658333333,7.475,6.383333333,6.2,5.241666667,4.658333333,4.7,5.9,6.8,6.9,6.241666667,5.566666667,5.058333333,4.375,4.175,3.966666667,4.116666667,5.191666667,5.683333333,5.283333333,4.458333333,4.125,4.066666667,4.591666667,6.883333333,11.725,13.73333333,13.31666667,11.60833333,9.966666667,8.158333333,6.85,5.808333333,5.016666667,4.408333333,3.9,13.01666667,8.333333333)
unemployment_cut <- unemployment[1:(44+1)]
une <- scale(unemployment_cut)

general <- read.csv(file = 'overall_hunting_data.csv')
general_nv <- subset(general, ST == 'NV' & Year > 1975)
res <- (general_nv$Resident.Licenses)


econ_data <- read.csv('economic_data.csv')
econ_data_sub <- subset(econ_data, Year > 1975)
PDI <- (econ_data_sub$Per.Capita.Personal.Disposable.Income/10000)
GAS <- econ_data_sub$Gas.Price..July.Oct.

mean_ratio <- 2.581635
sd_ratio <- 0.8894599

ratio <- PDI/GAS
############################################################################################################################
# Specify model
############################################################################################################################
code_gib <- nimbleCode( {
  
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(b0.gas + bt.gas[era[t]]*t, sd = sig.gas),0,)
    REL.COST[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  }
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  b0.pdi ~ dnorm(0, 0.001)
  bt.pdi ~ dnorm(0, 0.01)
  b0.gas ~ dnorm(1.5, 1)
  bt.gas[1] ~ dnorm(0, 0.01)
  bt.gas[2] ~ dnorm(0, 0.01)
  ar1 ~ dunif(-1,1)
  pdi.trend[1] <- b0.pdi + bt.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    mu.pdi[t] <- pdi.trend[t] + ar1 * ( PDI[t-1] - pdi.trend[t-1] )
    pdi.trend[t] <- b0.pdi + bt.pdi * t
  }
  
  for(l in 1:n.region){
    for(m in 1:n.species){
      sig[m,l] ~ dgamma(1,1)
      Delta[m,m,l] <- pow(Q[m,m,l], -0.5)
      Lambda[m,m,l] <- sig[m,l]
    }
    for (i in 2:n.species){
      for (j in 1:(i-1)){
        Lambda[i,j,l] <- 0
        Delta[i,j,l] <- 0
      }
    }
    Sigma[1:n.species,1:n.species,l] <- Lambda[1:n.species,1:n.species,l] %*% P[1:n.species,1:n.species,l] %*% Lambda[1:n.species,1:n.species,l]
    Q[1:n.species,1:n.species,l] ~ dinvwish(S = I[1:n.species,1:n.species,l], df = n.species + 1)
    P[1:n.species,1:n.species,l] <- Delta[1:n.species,1:n.species,l] %*% Q[1:n.species,1:n.species,l] %*% Delta[1:n.species,1:n.species,l]
    
    for (i in 1:n.species){
      for (j in 1:n.species){
        rho[i,j,l] <- Sigma[i,j,l]/sqrt(Sigma[i,i,l] * Sigma[j,j,l])
      }
    }
    for(k in 1:n.species){
      for(h in 1:n.year){
        pred1[k,l,h] <- inprod(beta.trend[k,l,1:K], ZZ[h,1:K,k,l])
        mu[k,l,h] <- lbo1[k,l] + inprod(beta.trend[k,l,1:K], ZZ[h,1:K,k,l]) + beta.general[k,l] * res[h] + beta.income[k,l] * REL.COST[h] + beta.drought2[k,l] * wpdsi[h,l] + beta.jobs[k,l] * une[h]
      }
    }
  }
  sig.une ~ dunif(0,5)
  for(h in 1:(n.year)){
    une[h] ~ dnorm(0, sd = sig.une)
  }
  for(l in 1:n.region){
    sig.wpdsi[l] ~ dunif(0,5)
    for(h in 1:n.year){
      wpdsi[h,l] ~ dnorm(0, sd = sig.wpdsi[l])
      hunt.eps[1:n.species,l,h] ~ dmnorm( mu[1:n.species,l,h], cov =  Sigma[1:n.species,1:n.species,l] )
    }
  }
  for (k in 1:n.species){
    for(l in 1:n.region){
      beta.drought2[k,l] ~ dnorm(mu.drought2[l], sd = sig.drought2[l])
      beta.jobs[k,l] ~ dnorm(mu.jobs[l], sd = sig.jobs[l])
      beta.general[k,l] ~ dnorm(mu.gen[l], sd  = sig.gen[l])
      beta.income[k,l] ~ dnorm(mu.incom[l], sd  = sig.incom[l])
      lbo1[k,l] ~ dnorm(5, sd = 1)
      log(H[k,1,l]) <- hunt.eps[k,l,1]
      z[k,1,l] ~  dpois(H[k,1,l])
      for(h in 2:(n.year)){
        log(H[k,h,l]) <- hunt.eps[k,l,h]
        z[k,h,l] ~  dpois(H[k,h,l])
      }
    }
  }
  for(l in 1:n.region){
    mu.incom[l] ~ dnorm(0, 0.01)
    mu.gen[l] ~ dnorm(0, 0.01)
    mu.drought2[l] ~ dnorm(0, 0.01)
    sig.gen[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.incom[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.drought2[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.jobs[l] ~ dnorm(0, 0.01)
    sig.jobs[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for (i in 1:n.species){
      sig.trend[i,l] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(k in 1:K){
        
        beta.trend[i,l,k] ~ dnorm(0, sd = sig.trend[i,l])
      }
    }
    for(h in 2:(n.year)){
      for (i in 1:n.species){
        lambda1[i,h-1,l] <- H[i,h,l]/H[i,h-1,l]
      }
    }
  }
  
})

# ###Matt's Changes
# code_gib <- nimbleCode( {
#     ### Predictors (Change in Hunters)
#     #Personal Disposable Income
#     b0.pdi ~ dnorm(0, 0.001) #intercept
#     bt.pdi ~ dnorm(0, 0.01) #year
#     sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
# 
#     ar1 ~ dunif(-1,1) #Change in PDI over time?
# 
#     pdi.trend[1] <- b0.pdi + bt.pdi * 1
#     mu.pdi[1] <- pdi.trend[1]
#     for(t in 2:n.year){
#       pdi.trend[t] <- b0.pdi + bt.pdi * t
#       mu.pdi[t] <- pdi.trend[t] + ar1 * (PDI[t-1] - pdi.trend[t-1])
#     } #t
# 
#     #Gas Prices
#     b0.gas ~ dnorm(1.5, 1)
#     sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
#     bt.gas[1] ~ dnorm(0, 0.01)
#     bt.gas[2] ~ dnorm(0, 0.01)
# 
#     #Relative Cost of Gas
#     for(t in 1:n.year){
#       PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
#       GAS[t] ~ T(dnorm(b0.gas + bt.gas[era[t]]*t, sd = sig.gas),0,)
#       REL.COST[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
#     } #t
# 
#     #Unemployment Rate
#     sig.une ~ dunif(0,5)
#     for(t in 1:(n.year)){
#       une[t] ~ dnorm(0, sd = sig.une)
#     } #t
# 
#     #Drought Index
#     for(r in 1:n.region){
#       sig.wpdsi[r] ~ dunif(0,5)
#       for(t in 1:n.year){
#         wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
#       } #t
#     } #r
# 
#     ### Priors for Regression Coefficients
#     for(r in 1:n.region){
#       mu.incom[r] ~ dnorm(0, 0.01)
#       sig.incom[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
#       mu.gen[r] ~ dnorm(0, 0.01)
#       sig.gen[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
#       mu.drought2[r] ~ dnorm(0, 0.01)
#       sig.drought2[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
#       mu.jobs[r] ~ dnorm(0, 0.01)
#       sig.jobs[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
# 
#       ### Coefficient Values
#       for(s in 1:n.species){
#         lbo1[s,r] ~ dnorm(5, sd = 1)
#         sig.trend[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
#         for(k in 1:K){
#           beta.trend[s,r,k] ~ dnorm(0, sd = sig.trend[s,r])
#         } #k
#         beta.drought2[s,r] ~ dnorm(mu.drought2[r], sd = sig.drought2[r])
#         beta.jobs[s,r] ~ dnorm(mu.jobs[r], sd = sig.jobs[r])
#         beta.general[s,r] ~ dnorm(mu.gen[r], sd  = sig.gen[r])
#         beta.income[s,r] ~ dnorm(mu.incom[r], sd  = sig.incom[r])
# 
#         ### Regression for Number Hunters
#         for(t in 1:n.year){
#           pred1[s,r,t] <- inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r])
#           mu[s,r,t] <- lbo1[s,r] + #Intercept
#             inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r]) + #Spline smoothing terms
#             beta.general[s,r] * res[t] + #Hunting Licences Sold
#             beta.income[s,r] * REL.COST[t] + #PDI/Gas Price
#             beta.drought2[s,r] * wpdsi[t,r] + #Drought
#             beta.jobs[s,r] * une[t] #Unemployment
#         } #t
#       } #s
# 
#       #Quantify annual change in number of hunters
#       for(t in 2:(n.year)){
#         for (s in 1:n.species){
#           lambda1[s,t-1,r] <- H[s,t,r]/H[s,t-1,r]
#         } #s
#       } #t
# 
# 
#       ### Define Variance-Covariance Matrix for Number of Hunters
#       Q[1:n.species,1:n.species,r] ~ dinvwish(S = I[1:n.species,1:n.species,r], df = n.species + 1)
#       for(s in 1:n.species){
#         sig[s,r] ~ dgamma(1,1)
#         Lambda[s,s,r] <- sig[s,r]
#         Delta[s,s,r] <- pow(Q[s,s,r], -0.5)
#       } #s
#       for (s1 in 2:n.species){
#         for (s2 in 1:(s1-1)){
#           Lambda[s1,s2,r] <- 0
#           Delta[s1,s2,r] <- 0
#         } #j
#       } #i
#       P[1:n.species,1:n.species,r] <- Delta[1:n.species,1:n.species,r] %*% Q[1:n.species,1:n.species,r] %*% Delta[1:n.species,1:n.species,r]
#       #Covariance matrix for multivariate normal dist. describing number of hunters by species
#       Sigma[1:n.species,1:n.species,r] <- Lambda[1:n.species,1:n.species,r] %*% P[1:n.species,1:n.species,r] %*% Lambda[1:n.species,1:n.species,r]
# 
# 
#       ### Process determining observed number of hunters
#       for(t in 1:n.year){
#         hunt.eps[1:n.species,r,t] ~ dmnorm(mu[1:n.species,r,t], cov =  Sigma[1:n.species,1:n.species,r])
#       } #t
# 
#       # Number of Hunters
#       for(s in 1:n.species){
#         for(t in 1:(n.year)){
#           log(H[s,t,r]) <- hunt.eps[s,r,t] #Log link to restrict H above 0
#           n.hunt[s,t,r] ~  dpois(H[s,t,r]) #Observed number of hunters species*year
#         } #t
#       } #s
# 
#       #Calculate Correlation Coefficient
#       for (s1 in 1:n.species){
#         for (s2 in 1:n.species){
#           #Covariance/((SD1*SD2)^.5)?
#           rho[s1,s2,r] <- Sigma[s1,s2,r]/sqrt(Sigma[s1,s1,r] * Sigma[s2,s2,r])
#         } #s2
#       } #s1
# 
#     } #r
# 
#   })

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

library(splines)
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


library(LaplacesDemon)

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

Sigma = nearPD(Sigma, corr = FALSE,doSym = TRUE)


constants_gib <- list(      n.species = 7, K= 12,
                            n.region = 2,
                            n.year = cut+4,
                            mu.hunt = rep(0, 7),
                            era = c(rep(1,19),rep(2, 27)),
                            mean.H = apply(hunters, c(1,3), mean, na.rm = TRUE),
                            sd.H = apply(hunters, c(1,3), sd, na.rm = TRUE),
                            I = abind(I,I,along = 3))

data_gib <- list(une = c(une,NA), PDI = PDI, GAS = GAS,
                 z = abind(hunters,array(NA, dim = c(7,4,2)) ,along = 2), 
                 ZZ = ZZ1, res= scale(res)[,1],
                 wpdsi = data.matrix(abind(wpdsi,matrix(NA,1,2), along = 1)))

# Initial values
Hi <- hunters + 2500
Hi[,-1,] <- NA
Hi[7,10,] <- rpois(2,colMeans(hunters[7,-10,]))
Hi <- abind(Hi[,1:cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi <- abind(hunters,hunters[,cut,],hunters[,cut,],hunters[,cut,],hunters[,cut,], along = 2)
zi[7,10,] <- Hi[7,10,]  

rho.hunt.init <- diag(7)
PDI.inits <- ifelse(is.na(PDI) == TRUE, mean(PDI, na.rm = TRUE), PDI)
GAS.inits <- ifelse(is.na(GAS) == TRUE, 0.9, GAS)
# Initial values
initsFunction <- function() list( 
  GAS = GAS.inits, PDI = PDI.inits, sig.pdi = 1, sig.gas = 1,
  beta.drought2 = matrix(0, 7, 2), mu.drought2 = c(0,0), sig.drought2 = c(1,1),
  beta.income = matrix(0, 7, 2), mu.incom = c(0,0), sig.incom= c(1,1),
  beta.jobs = matrix(0, 7, 2), mu.jobs = c(0,0), sig.jobs = c(1,1), sig.wpdsi = c(1,1), sig.une = 1,
  
  b0.pdi = 2.9, bt.pdi = 0, b0.gas = 0.9,  ar1 = 0, bt.gas = c(0,0),
  Q = abind(Q,Q,along = 3),
  Sigma = abind(Sigma,Sigma,along = 3),
  P = abind(P,P,along = 3),
  Lambda = abind(Lambda,Lambda,along = 3),
  Delta = abind(Delta,Delta,along = 3),
  rho = abind(diag(n.species),diag(n.species),along = 3),
  
  H = Hi, x = zi,
  sig.trend = matrix(1, ncol = 2, nrow = 7),
  beta.general = matrix(0, ncol = 2, nrow = 7), 
  mu.gen = c(0,0), sig.gen = c(1,1),
  lbo1 = matrix(0, ncol = 2, nrow = 7),
  
  lambda1 = array(1, dim = c(7,cut+3,2) ),
  log.r1  = array(0, dim = c(7,cut+3,2) ),
  hunt.eps = array(rnorm(7*2*cut+4,0,0.1),dim = c(7,2,cut+4)),
  
  wpdsi = matrix(0,46,2), WPDSI = matrix(0,46,2),
  une = rep(0,46), UNE = rep(0,46)
  
  
)

# Parameters monitored
pars1 <- c('b0')

start_time <- Sys.time()

inits_gib <- initsFunction()

library(parallel)
library(coda)

nc <- 3    # number of chains
cl<-makeCluster(nc,timeout=5184000)
clusterExport(cl, c("code_gib", "inits_gib", "data_gib", "constants_gib", "pars1"))
for (j in seq_along(cl)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl[j], "inits")
}

out_gib <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  model_test <- nimbleModel( code = code_gib, constants = constants_gib,  data =  data_gib, inits = inits_gib )
  
  model_test$simulate(c('sig', 'pred1', 'mu', 'une', 'wpdsi',  'beta.trend'))
  model_test$initializeInfo()
  model_test$calculate()
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  c('lbo1','beta.drought2','beta.jobs','wpdsi', 'une','GAS','PDI', 'REL.COST',
                                                           'hunt.eps', 'H','mu','rho','pred1','ar1', 'beta.income',
                                                           'lambda1', 'beta.trend', 'beta.general'))
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 4000, niter = 6000, thin = 1, thin2 = 1)
  
  return(samplesList)
})

end_time <- Sys.time()
end_time - start_time

stopCluster(cl)

samples2 <- list(chain1 =  out_gib[[1]]$samples2, 
                 chain2 =  out_gib[[2]]$samples2, 
                 chain3 =  out_gib[[3]]$samples2)

samples1    <- list(chain1 =  out_gib[[1]]$samples, 
                    chain2 =  out_gib[[2]]$samples, 
                    chain3 =  out_gib[[3]]$samples)

library(coda)
library(MCMCvis)

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

test1 <- MCMCsummary(mcmcList2, 'pred1')
test1.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)],  region = c('east','west'), year = 1976:2021 ), test1)

fig1 <- ggplot(data = test1.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Numbers', x = 'Year') + facet_grid(region ~ species)
fig1

test3 <- MCMCsummary(mcmcList2, 'pred2')
test3.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)],  region = c('east','west'), year = 1976:2021 ), test3)

fig3 <- ggplot(data = test3.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Numbers', x = 'Year') + facet_grid(region ~ species) +  scale_fill_flat_d() + theme_modern()
fig3



hunter_prime <- abind(hunters, array(NA, dim = c(7,4,2)), along = 2)
test2    <- MCMCsummary(mcmcList2, 'H')
test2.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], year = 1976:2021, region = c('east','west') ), test2, hunters = c(hunter_prime))

fig2 <- ggplot(data = test2.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_point(aes(x =year, y = hunters )) +
  geom_hline(yintercept = 0) + labs(y = 'Hunter Numbers', x = 'Year') + facet_wrap(region ~ species, ncol  = 7, scales = 'free') +  scale_fill_flat_d() + theme_modern()
fig2

############################################################################################################################
# Specify model
############################################################################################################################

hunters.prime   <- MCMCpstr(mcmcList2, 'H')$H
nseg <- 10
BM <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z  <- array(NA, dim = c(cut+4,nseg+2,7,2))
D <- diff(diag(ncol(BM[,,1,1])), diff = 1)
Q <- t(D) %*% solve(D %*% t(D))

for(i in 1:7){
  for(j in 1:2){
    BM[,,i,j] <- bs_bbase(hunters.prime[i,,j], nseg = 10)
    Z[,,i,j] <-  BM[,,i,j]%*% Q
  }
}

ZZ <- Z
ZZ[is.na(ZZ)] <- 0

code <- nimbleCode( {
  for(l in 1:n.region){
    for(m in 1:n.species){
      sig[m,l] ~ dgamma(1,1)
      Delta[m,m,l] <- pow(Q[m,m,l], -0.5)
      Lambda[m,m,l] <- sig[m,l]
      
      sig2[m,l] ~ dgamma(1,1)
      Delta2[m,m,l] <- pow(Q2[m,m,l], -0.5)
      Lambda2[m,m,l] <- sig2[m,l]
    }
    for (i in 2:n.species){
      for (j in 1:(i-1)){
        Lambda[i,j,l] <- 0
        Delta[i,j,l] <- 0
        
        Lambda2[i,j,l] <- 0
        Delta2[i,j,l] <- 0
      }
    }  
    Sigma[1:n.species,1:n.species,l] <- Lambda[1:n.species,1:n.species,l] %*% P[1:n.species,1:n.species,l] %*% Lambda[1:n.species,1:n.species,l]  
    Q[1:n.species,1:n.species,l] ~ dinvwish(S = I[1:n.species,1:n.species,l], df = n.species + 1)
    P[1:n.species,1:n.species,l] <- Delta[1:n.species,1:n.species,l] %*% Q[1:n.species,1:n.species,l] %*% Delta[1:n.species,1:n.species,l]
    
    Sigma2[1:n.species,1:n.species,l] <- Lambda2[1:n.species,1:n.species,l] %*% P2[1:n.species,1:n.species,l] %*% Lambda2[1:n.species,1:n.species,l]  
    Q2[1:n.species,1:n.species,l] ~ dinvwish(S = I2[1:n.species,1:n.species,l], df = n.species + 1)
    P2[1:n.species,1:n.species,l] <- Delta2[1:n.species,1:n.species,l] %*% Q2[1:n.species,1:n.species,l] %*% Delta2[1:n.species,1:n.species,l]  
    
    for (i in 1:n.species){
      for (j in 1:n.species){
        rho[i,j,l] <- Sigma[i,j,l]/sqrt(Sigma[i,i,l] * Sigma[j,j,l])   
        rho2[i,j,l] <- Sigma2[i,j,l]/sqrt(Sigma2[i,i,l] * Sigma2[j,j,l])   
      }
    }
    
    sig.wpdsi[l] ~ dunif(0,5)
    sig.pdsi[l] ~ dunif(0,5)
    for(h in 1:n.year){
      wpdsi[h,l] ~ dnorm(0, sd = sig.wpdsi[l])
      pdsi[h,l] ~ dnorm(0, sd = sig.pdsi[l])
    }        
    
    for(k in 1:n.species){
      for(h in 1:n.year){
        pred1[k,l,h] <- inprod(beta.trend[k,l,1:K], ZZ[h,1:K,k,l])
        mu[k,l,h] <- lbo1[k,l] + inprod(beta.trend[k,l,1:K], ZZ[h,1:K,k,l]) + beta.drought2[k,l] * wpdsi[h,l] + beta.jobs[k,l] * une[h]
      }
    }
  }
  sig.une~ dunif(0,5)
  for(h in 1:(n.year)){
    une[h] ~ dnorm(0, sd = sig.une)
    for(l in 1:n.region){
      hunt.eps[1:n.species,l,h] ~ dmnorm( mu[1:n.species,l,h], cov =  Sigma[1:n.species,1:n.species,l] )
    }
  }
  for (k in 1:n.species){
    for(l in 1:n.region){
      beta.drought2[k,l] ~ dnorm(mu.drought2[l], sd = sig.drought2[l]) 
      beta.jobs[k,l] ~ dnorm(mu.jobs[l], sd = sig.jobs[l]) 
      lbo1[k,l] ~ dnorm(5, sd = 1)
      lbo[k,l] ~ dnorm(0, sd = 1)
      log(H[k,1,l]) <- hunt.eps[k,l,1]
      z[k,1,l] ~  dpois(H[k,1,l])
      for(h in 2:(n.year)){
        log(H[k,h,l]) <- hunt.eps[k,l,h]
        z[k,h,l] ~  dpois(H[k,h,l])
      }
    }
  }
  for(l in 1:n.region){
    mu.drought[l] ~ dnorm(0, 0.01)
    sig.drought[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.drought2[l] ~ dnorm(0, 0.01)
    sig.drought2[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.jobs[l] ~ dnorm(0, 0.01)
    sig.jobs[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for (i in 1:n.species){
      beta.drought[i,l] ~ dnorm(mu.drought[l], sd = sig.drought[l])
      sig.pressure[i,l] ~ T(dt(0, pow(2.5,-2), 1),0,)
      sig.trend[i,l] ~ T(dt(0, pow(2.5,-2), 1),0,)
      for(k in 1:K){
        beta.pressure[i,l,k] ~ dnorm(0, sd = sig.pressure[i,l]) 
        beta.trend[i,l,k] ~ dnorm(0, sd = sig.trend[i,l])
      }
    } 
    for (i in 1:n.species){
      N[i,1,l] ~ dpois(y[i,1,l])
      bph[i,1,l] <-  N[i,1,l]/H[i,1,l]
      for(h in 2:(n.year)){
        mu2[i,h-1,l] <- lbo[i,l] + inprod(beta.pressure[i,l,1:K], Z[h-1,1:K,i,l]) + beta.drought[i,l] * pdsi[h-1,l]
      }
    }
    for(h in 2:(n.year)){
      log.r[1:n.species,h-1,l]  ~ dmnorm( mu2[1:n.species,h-1,l], cov =  Sigma2[1:n.species,1:n.species,l])
    }
    for(h in 2:(n.year)){
      for (i in 1:n.species){
        pred[i,h-1,l] <- inprod(beta.pressure[i,l,1:K], Z[h-1,1:K,i,l]) 
        lambda[i,h-1,l] <- exp(log.r[i,h-1,l])
        lambda1[i,h-1,l] <- H[i,h,l]/H[i,h-1,l]
        N[i,h,l] <- lambda[i,h-1,l] * N[i,h-1,l]
        y[i,h,l] ~  dpois(N[i,h,l]) 
        bph[i,h,l] <-  N[i,h,l]/H[i,h,l]
      }
    }
  }
  
  for(l in 1:n.site){
    sigma.chuk[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    lbo2[l] ~ dlogis(0,1)
    theta2[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    C[l,1] ~ dpois(x[l,1])
    mod[l] ~ dlogis(0,1)
    for(h in 2:n.yr){
      log.r2[l,h-1] <- lbo2[l] + chuk.eps[l,h-1]
      lambda2[l,h-1] <- exp(log.r2[l,h-1])
      C[l,h] <- lambda2[l,h-1] * C[l,h-1]
      rate2[l,h] <- theta2[l]/(theta2[l] + C[l,h])
      x[l,h] ~   dnegbin(rate2[l,h], theta2[l])  
    }
    for(h in 1:(n.yr-1)){
      chuk.eps[l,h]  ~ dnorm(mod[l] * log.r[3,h+14,1], sd = sigma.chuk[l])
    }
  }
  
  
})


Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

library(splines)
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
BM <- array(NA, dim = c(cut+4,nseg+3,7,2))
Z  <- array(NA, dim = c(cut+4,nseg+2,7,2))
D <- diff(diag(ncol(BM[,,1,1])), diff = 1)
Q <- t(D) %*% solve(D %*% t(D))

hunter.prime <- abind(hunters, array(NA, dim = c(7,4,2)),along = 2)

for(i in 1:7){
  for(j in 1:2){
    BM[,,i,j] <- bs_bbase(hunter.prime[i,1:(cut+4),j], nseg = 10)
    Z[,,i,j] <-  BM[,,i,j]%*% Q
  }
}

ZZ <- Z
ZZ[is.na(ZZ)] <- 0

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

constants <- list(      n.species = 7, K= 12,
                        n.region = 2,
                        n.year = cut+4,
                        mu.hunt = rep(0, 7),
                        mu.bird = rep(0, 7),
                        n.site = 13, diff = 14,
                        n.yr = ncol(chukar),
                        mean.H = apply(hunter.prime, c(1,3), mean, na.rm = TRUE),
                        sd.H = apply(hunter.prime, c(1,3), sd, na.rm = TRUE),
                        I = abind(I,I,along = 3), I2 =abind(I2,I2,along = 3))

data <- list(y = abind(upland,array(NA, dim = c(7,4,2)) ,along = 2),  une = c(une,NA),
             z = abind(hunters,array(NA, dim = c(7,4,2)) ,along = 2), 
             Z = ZZ, ZZ = ZZ1, 
             pdsi = data.matrix(abind(pdsi,matrix(NA,2,2), along = 1)),
             wpdsi = data.matrix(abind(wpdsi,matrix(NA,1,2), along = 1)), 
             x =  data.matrix(chukar))

# Initial values
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

# Initial values
initsFunction <- function() list( 
  phi = matrix(rbeta(14,1,1),7,2),
  beta.drought = matrix(0, 7, 2), mu.drought = c(0,0), sig.drought = c(1,1),
  beta.drought2 = matrix(0, 7, 2), mu.drought2 = c(0,0), sig.drought2 = c(1,1),
  beta.jobs = matrix(0, 7, 2), mu.jobs = c(0,0), sig.jobs = c(1,1),
  
  C = chukar_na + 50,
  x =  chukar_na,
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
  N = Ni,  H = Hi, z= zi, y = yi,
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

# Parameters monitored
pars1 <- c('b0')

start_time <- Sys.time()

inits <- initsFunction()


library(parallel)
library(coda)
set.seed(2)
nc <- 3    # number of chains
cl<-makeCluster(nc,timeout=5184000)
clusterExport(cl, c("code", "inits", "data", "constants", "pars1"))
for (j in seq_along(cl)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl[j], "inits")
}

out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  model_test <- nimbleModel( code = code, constants = constants,  data =  data, inits = inits )
  
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  c('lbo1','lbo','N', 'lambda','pred','beta.drought2','beta.drought','beta.jobs',
                                                           'hunt.eps', 'H','bph','mu','mu2','pred1','rho','rho2',
                                                           'lambda1', 'beta.pressure','beta.trend',
                                                           'lbo2', 'chuk.eps', 'lambda2',
                                                           'C', 'lambda2','theta2','log.r2', 'sigma.chuk'))
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  samplesList <- runMCMC(Cmcmc,nburnin = 250000, niter = 500000, thin = 10, thin2 = 10)
  
  return(samplesList)
})

end_time <- Sys.time()
end_time - start_time

stopCluster(cl)

samples2 <- list(chain1 =  out[[1]]$samples2, 
                 chain2 =  out[[2]]$samples2, 
                 chain3 =  out[[3]]$samples2)

samples1    <- list(chain1 =  out[[1]]$samples, 
                    chain2 =  out[[2]]$samples, 
                    chain3 =  out[[3]]$samples)

library(coda)
library(MCMCvis)

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))



mcmcList1 <- as.mcmc.list(lapply(list(samplesList $samples), mcmc))
mcmcList2 <- as.mcmc.list(lapply(list(samplesList $samples2), mcmc))

files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'model_output.rdata')

rho1 <- MCMCpstr(mcmcList2,'rho')
rho2 <- MCMCpstr(mcmcList2,'rho2')

rho1a <- rho1$rho[,,1]
rho1b <- rho1$rho[,,2]
rownames(rho1a) <- colnames(rho1a) <- species[-c(4,8)]
rownames(rho1b) <- colnames(rho1b) <- species[-c(4,8)]

rho1a[lower.tri(rho1a)] <- rho1b[lower.tri(rho1b)]

rho2a <- rho2$rho2[,,1]
rho2b <- rho2$rho2[,,2]
rownames(rho2a) <- colnames(rho2a) <- species[-c(4,8)]
rownames(rho2b) <- colnames(rho2b) <- species[-c(4,8)]

rho2a[lower.tri(rho2a)] <- rho2b[lower.tri(rho2b)]

diag(rho1a) <-NA
diag(rho2a) <-NA

library(ggcorrplot)

ggcorrplot(rho1a, outline.col = "white",  show.diag = FALSE, lab= TRUE)
ggcorrplot(rho2a, outline.col = "white",  show.diag = FALSE, lab= TRUE)

summary_one <- MCMCsummary(mcmcList1)
summary_two <- MCMCsummary(mcmcList2)

test.bph <- MCMCsummary(mcmcList2, 'bph')
test.H   <- MCMCsummary(mcmcList2, 'H')
test.N   <- MCMCsummary(mcmcList2, 'N')

test.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                        year = 1976:2021, 
                                        region = c('west','east') ),
                            test.H, hunters = c(abind(hunters, array(NA, dim = c(7,4,2)),along = 2)))

test.dfn <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                         year = 1976:2021, 
                                         region = c('west','east') ),
                             test.N, upland = c(abind(upland, array(NA, dim = c(7,4,2)),along = 2)))

fig_upland <-   ggplot(data = test.dfn, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`),  size = .75) +
  geom_pointrange2(data = subset(test.dfn,year>2017), aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), color = 'red',  size = .75) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1000, color = NA) +  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested', x = 'Year') +
  scale_fill_flat_d() + theme_modern()
fig_upland

fig_hunters <-  ggplot(data = test.df, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`),  size = .75) +
  geom_pointrange2(data = subset(test.df,year>2017), aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), color = 'red',  size = .75) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 0) + labs(y = 'Predicted number of hunters', x = 'Year') +
  scale_fill_flat_d() + theme_modern()
fig_hunters

ggsave(fig_hunters, file = 'fig_hunters.tiff', height = 6, width = 15, dpi = 320)

ggplot(data = test.df, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), fill = 'black', shape = 21, size = .75, position = position_dodge(1)) +
  geom_point(aes(x = year, y = hunters), color = 'red') +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1) +
  scale_fill_flat_d() + theme_modern()

test <- MCMCsummary(mcmcList2, 'pred')
test.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], year = 1977:2021, region = c('east','west') ),
                            test, hunters =  c(abind(hunters, array(NA, dim = c(7,3,2)),along = 2)))

fig1 <- ggplot(data = test.df, aes(x = hunters, y = mean, group = region)) +
  geom_pointrange2(aes(x =hunters, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Pressure', x = 'Year') + facet_wrap(region ~ species, scales = 'free', ncol = 7) +  scale_fill_flat_d() + theme_modern(legend.position = 'none')
fig1

test1 <- MCMCsummary(mcmcList2, 'pred1')
test1.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)],  region = c('east','west'), year = 1976:2021 ), test1)

fig2 <- ggplot(data = test1.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Numbers', x = 'Year') + facet_grid(region ~ species) +  scale_fill_flat_d() + theme_modern()
fig2

test.NH <- cbind.data.frame(expand.grid(species = species[-c(4,8)], year = 1976:2021, region = c('west','east')), test.bph)

fig1 <- ggplot(data = subset(test.NH, year == 2017), aes(x = species, y = mean, group = region)) +
  geom_pointrange2(aes(x =species, y = mean, ymin = `2.5%`,ymax = `97.5%`, shape = region),color = 'black', size = .75, position = position_dodge(1)) +
  # geom_point(aes(x = species, y = bph,shape = region),size = 2, color = 'red', position = position_dodge(1)) +
  geom_text(data = subset(test.dfn, year == 2017), aes(x =species, y = -1, group = region, label = round(mean,0)),position = position_dodge(1))+
  geom_text(data = subset(test.dfn, year == 2017), aes(x =species, y = -3, group = region, label = round(upland,0)),position = position_dodge(1), color = 'red')+
  geom_hline(yintercept = 0) + labs(y = 'Forecasted number of individuals harvested per hunter in 2017', x = 'Species') +
  theme_modern()
fig1
ggsave(fig1, file = 'prelim_model.tiff', dpi = 320, height = 6, width = 10)

ggplot(data = test.NH, aes(x = year, y = mean, fill = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, shape = region),color = 'black', size = .75, position = position_dodge(1)) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1) +
  scale_fill_flat_d() + theme_modern()

ggplot(data = test.NH, aes(x = year, y = mean, fill = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`), shape = 21, size = .75, position = position_dodge(1)) +
  facet_wrap(region ~ species, ncol = 7, scales = 'free_y') +
  geom_hline(yintercept = 0) + labs(y = 'Estimated number of individuals harvested per hunter', x = 'Species') +
  scale_fill_flat_d() +  scale_color_metro_d() + theme_modern()

library(see)
ggplot(data = test.df, aes(x = hunters, y = mean, fill = region)) +
  geom_pointrange2(aes(x = hunters, y = mean, ymin = `2.5%`,ymax = `97.5%`), linetype = 'dashed', shape = 21, size = .66, position = position_dodge(1)) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  labs(x = 'Number of self-reported active hunters in year t', y = 'Observed impact on the number of individuals reportd as harvested in year t+1') +
  geom_hline(yintercept = 0) +
  scale_fill_flat_d() + theme_modern()

ggplot(data = test.df, aes(x = year, y = mean, fill = region)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), linetype = 'dashed', shape = 21, size = .66, position = position_dodge(1)) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  labs(x = 'Year', y = 'Observed impact of active hunters in year t on the number\n of individuals reported as harvested in year t+1') +
  geom_hline(yintercept = 0) +
  scale_fill_flat_d() + theme_modern()

plot(hunters[1,1:41,1], test_matrix[1,,1])


N_bird   <- MCMCsummary(mcmcList2, 'N')
N_hunt   <- MCMCsummary(mcmcList2, 'H')
N_survey <- MCMCsummary(mcmcList2, 'C')

N_bird$Species <- rep(species[-c(4,8)], times = 46 * 2)
N_bird$Year    <- rep(1976:2021, each = 7, times = 2)
N_bird$Region  <- rep(c("West",'East'), each = 7 *  46)
N_bird$Type <- rep('Bird')

N_hunt$Species <- rep(species[-c(4,8)], times = 46 * 2)
N_hunt$Year    <- rep(1976:2021, each = 7, times = 2)
N_hunt$Region  <- rep(c("West",'East'), each = 7 *  46)
N_hunt$Type <- rep('Hunter')


hunting <- rbind(N_bird, N_hunt)

SGN <- subset(hunting, grepl('SAGR', Species))

dodge = position_dodge(.5)
ggplot(SGN, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Sage-grouse Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

CHN <- subset(hunting, grepl('CHUK', Species))

dodge = position_dodge(.5)
ggplot(CHN, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Chukar Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

BLG <- subset(hunting, grepl('BLGR', Species))

dodge = position_dodge(.5)
ggplot(BLG, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Blue-grouse Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

HUP <- subset(hunting, grepl('HUPA', Species))

dodge = position_dodge(.5)
ggplot(HUP, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Hungarian Partridge Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

CAQ <- subset(hunting, grepl('CAQU', Species))

dodge = position_dodge(.5)
ggplot(CAQ, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and California Quail Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

PHE <- subset(hunting, grepl('PHEA', Species))

dodge = position_dodge(.5)
ggplot(PHE, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region)) + #, formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Ring-necked Pheasant Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')


RAB <- subset(hunting, grepl('RABB', Species))

dodge = position_dodge(.5)
ggplot(RAB, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region)) + #, formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Rabbit Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')


lambda_bird <- MCMCsummary(mcmcList2, 'lambda')
lambda_hunt <- MCMCsummary(mcmcList2, 'lambda1')
lambda_survey <- MCMCsummary(mcmcList2, 'lambda2')


species_mod <- species[-c(4,8)]

lambda_bird$Species <- rep(species_mod, times = 41 * 2)
lambda_bird$Year   <- rep(1977:2017, each = 7 * 2)
lambda_bird$Region <- rep(c("West",'East'), each = 7, times = 41)
lambda_bird$Type <- rep('Bird')

ggplot(lambda_bird, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, fill = Region), shape  =21, color = 'black') +
  scale_fill_flat_d() +
  facet_wrap(~Species, ncol = 2, scales = 'free_y') + theme_bw() + theme_modern()

lambda_hunt$Species <- rep(species_mod, times = 41 * 2)
lambda_hunt$Year   <- rep(1977:2017, each = 7 * 2)
lambda_hunt$Region <- rep(c("West",'East'), each = 7, times = 41)
lambda_hunt$Type <- rep('Hunter')

ggplot(lambda_hunt, aes(x = Year, y = log(`50%`), color = Region)) +
  geom_pointrange(aes(x = Year, ymin = log(`2.5%`), ymax = log(`97.5%`), fill = Region), shape  =21, color = 'black') +
  scale_fill_flat_d() +
  facet_wrap(~Species, ncol = 3, scales = 'free_y') + theme_bw() + theme_modern()



lambdas <- rbind(lambda_bird, lambda_hunt)

SG <- subset(lambdas, grepl('SAGR', Species))

dodge = position_dodge(.5)
ggplot(SG, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()


CH <- subset(lambdas, grepl('CHUK', Species))
ggplot(CH, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','black')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()


CQ <- subset(lambdas, grepl('CAQU', Species))
ggplot(CQ, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','black')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()


HP <- subset(lambdas, grepl('HUPA', Species))
ggplot(HP, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','black')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()
