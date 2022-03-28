### Chukar Site Abundance
# Breaking down model to identify specific relationships

set.seed(2)
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "coda", "MCMCvis", "parallel"),
       require, character.only = T)

### Model
code <- nimbleCode( {
  ### Chukar Site Abundance
  for(r in 1:n.reg){
    sig.rabbit[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    beta.rabbit[r] ~ dnorm(0, sd = sig.rabbit[r])
  }
  
  for(p in 1:n.site){
    #Regression priors
    lbo2[p] ~ dlogis(0,1) #Intercept
    sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mod[p] ~ dlogis(0,1)
    
    #NB priors
    theta2[p] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
    C[p,1] ~ dpois(x[p,1]) #Equivalent of Poisson lambda
    
    for(t in 1:(n.yr-1)){
      #Regression
      chuk.eps[p,t]  ~ dnorm(mod[p], sd = sigma.chuk[p]) #Coeffient for annual change
      log.r2[p,t] <- lbo2[p] + beta.rabbit[site.reg[p]]*rabbit.h[t, site.reg[p]] + chuk.eps[p,t] #Unlinked change in abundance
      
      #NB for Count
      lambda2[p,t] <- exp(log.r2[p,t]) #Change in abundance
      C[p,t+1] <- lambda2[p,t] * C[p,t] #Equivalent of Poisson lambda
      rate2[p,t+1] <- theta2[p]/(theta2[p] + C[p,t+1]) #NB success parameter
      n.chuk[p,t+1] ~ dnegbin(rate2[p,t+1], theta2[p]) #obs. # of chukars follow neg-bin  
    } #t
  } #p  
})


### Data
## Site Counts
survey_data <- read.csv('./Data/Chukar_Surveys_data.csv') #Number of surveys for each population per year
df <- gather(survey_data, "Population", "Count", 2:14) # Where chukar surveys occurred
site_list <- unique(df$Population) # Years to assess
year_list <- rep(1972:2017) #Create new df for all population*year combinations
df_full <- data.frame("Population" = rep(site_list, times = length(year_list)), "Year" = rep(year_list, each = length(site_list))) #Join associated number of surveys at each location in each year
df_full <- left_join(df_full, df)
df_ch <- dcast(df_full, Population ~ Year)
chukar <- df_ch[,-c(1:19)]

## Site Region Association
site_region <- read.csv('./Data/Chukar_Surveys_locations.csv') %>%
  select(Population = Survey.Location, Region = NDOWREGION) %>%
  mutate(Index = ifelse(Region == "Eastern", 1, 2))
site_region.df <- merge(df_ch, site_region, by = "Population") %>% select(Population, Region, Index)
site_region <- site_region.df[,3]

## Rabbit Harvest
rabbits <- as.matrix(read.csv('./Data/all_species_harvest_data.csv') %>%
                       filter(Species == "RABBIT", Section !="Southern") %>%
                       select(Region = Section, Year, Rabbits = Animals) %>%
                       pivot_wider(values_from = "Rabbits", names_from = "Region") %>%
                       select(Eastern, Western))


## Wrap in List
data <- list(n.chuk =  data.matrix(chukar),
             rabbit.h = rabbits
             )

### Initial Values
chukar_na <- chukar 
chukar_na <- ifelse(is.na(chukar == TRUE), 1, 0)

initsFunction <- function() list( 
  C = chukar_na + 50,
  n.chuk =  chukar_na,
  sigma.chuk = rep(1,13),
  theta2 = rep(1,13),
  lbo2 =  rep(0,13)
)

inits <- initsFunction()

### Contansts
constants <- list(
                  n.site = 13,
                  n.yr = ncol(chukar),
                  n.reg = 2,
                  site.reg = site_region
                  )

### Monitors
pars1 <- c(
           'lbo2', 
           'chuk.eps',
           'C', 
           'lambda2',
           'theta2',
           'log.r2', 
           'sigma.chuk',
           'sig.rabbit',
           'beta.rabbit'
           )


### MCMC
#Check Model
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits)

model_test$simulate(pars1)
model_test$initializeInfo()
model_test$calculate()

#Parallel Processing Setup
nc <- detectCores()/2    # number of chains
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
  model_test <- nimbleModel( code = code, 
                             constants = constants,  
                             data =  data, 
                             inits = inits )
  
  mcmcConf <-  configureMCMC( model_test,   monitors2 =  pars1)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  
  # samplesList <- runMCMC(Cmcmc,nburnin = 250000, niter = 500000, thin = 10, thin2 = 10)
  samplesList <- runMCMC(Cmcmc,nburnin = 2500, niter = 5000, thin = 10, thin2 = 10)
  
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



mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))



mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

files <- list(mcmcList1,mcmcList2,code)
save(files, file = 'ChukarSEM_model_output.rdata')