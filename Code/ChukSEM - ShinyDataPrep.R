lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr", "tidyr", "mvtnorm", "R.utils"), require, character.only = T)


#Load model outputs to get beta and sigma estimates
load(file = "./output/NDOW_Upland_SEM_output.rdata")
beta.output <- files[[2]]
rho.output <- files[[1]]
data.input <- files[[4]]

appobject <- list()

appobject$N.data <- data.input$n.harv
appobject$H.data <- data.input$n.hunt

### Hunter Effort Covariance Matrix
rho.harv.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 5, ncol = 5))
appobject$rho.harv <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.harv.est)){
  appobject$rho.harv[[rho.harv.est$Region[i]]][rho.harv.est$Species1[i], rho.harv.est$Species2[i]] <- rho.harv.est$Estimate[i]
}

### Hunter Effort Covariance Matrix
rho.hunt.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 5, ncol = 5))
appobject$rho.hunt <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.hunt.est)){
  appobject$rho.hunt[[rho.hunt.est$Region[i]]][rho.hunt.est$Species1[i], rho.hunt.est$Species2[i]] <- rho.hunt.est$Estimate[i]
}

### Observed Hunter and Harvest Data For Graphs
appobject$H <- H.input <- data.input$n.hunt
appobject$N <- N.input <- data.input$n.harv


### Hunter Effort Regression Coefficients
Hunt.reg.coef1 <- MCMCsummary(beta.output, 'alpha.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID))%>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "b.0") %>%
  dplyr::select(Species, Region, A.hunt = mean, A.hunt.sd = sd) 
appobject$H.reg <- Hunt.reg.coef <- MCMCsummary(beta.output, 'beta.econ.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.ECON") %>%
  dplyr::select(Species, B.econ.hunt = mean, B.econ.hunt.sd = sd) %>%
  merge(., expand.grid(Species = 1:5, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Hunt.reg.coef1, ., by = c("Species", "Region"), all =T) %>%
  group_split(Region, .keep = F)

latent.trend <-list()
latent.trend$mean <- MCMCsummary(beta.output, 'latent.trend') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "Latent") %>%
  dplyr::select(Species, Region, Year, mean) %>%
  pivot_wider(names_from = Year, names_prefix = "Year", values_from = c(mean)) %>%
  dplyr::select(-Species) %>%
  group_split(Region, .keep = F)
latent.trend$sd <- MCMCsummary(beta.output, 'latent.trend') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "Latent") %>%
  dplyr::select(Species, Region, Year, sd) %>%
  pivot_wider(names_from = Year, names_prefix = "Year", values_from = c(sd))  %>%
  dplyr::select(-Species) %>%
  group_split(Region, .keep = F)
appobject$latent.trend <- latent.trend


### Total Harvest Regression Coefficients
Harv.reg.coef1 <- MCMCsummary(beta.output, 'alpha.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "alpha.harv") %>%
  dplyr::select(Species, Region, A.harv = mean, A.harv.sd = sd) 
  
Harv.reg.coef2 <- MCMCsummary(beta.output, 'beta.wintsev.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.WINT") %>%
  dplyr::select(Species, B.ws.harv = mean, B.ws.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:5, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef1, ., by = c("Species", "Region"), all =T)

Harv.reg.coef3 <- MCMCsummary(beta.output, 'beta.pdsi.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.PDSI") %>%
  dplyr::select(Species, B.pdsi.harv = mean, B.pdsi.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:5, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef2, ., by = c("Species", "Region"), all =T)

Harv.reg.coef4 <- MCMCsummary(beta.output, 'beta.bbs.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.BBS") %>%
  dplyr::select(Species, B.bbs.harv = mean, B.bbs.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:5, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef3, ., by = c("Species", "Region"), all =T)
  
appobject$N.reg <- Harv.reg.coef <- MCMCsummary(beta.output, 'beta.hunter.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(Species, Region, B.hunter = mean, B.hunter.sd = sd) %>%
  merge(Harv.reg.coef4, ., by = c("Species", "Region"), all =T) %>%
  group_split(Region, .keep = F)


save(appobject, file = "./NDOWSEM/www/NDOW_SEM_app_objects.R")

# #### For Integration into Shiny
# #These will represent the slider inputs. Will need to be reactive
# require(nimble)
# require(mvtnorm)
# 
# econ.var <- runif(1,-2,2)
# winter.var <- runif(1,-2,2)
# drought.var <- runif(1,-2,2)
# predator.var <- runif(1,-2,2)
# start.year.f <- 2015 - 1975
# y.b <- start.year.f-1 #years of data before forecast
# last.year.f <- 2017 - 1975
# n.year.f <- length(start.year.f:last.year.f) #If you want to show more than next year, use this as another dimension in array
# n.samples <- 1000
# 
# ### Hunter Effort
# H.reg <- appobject$H.reg
# N.reg <- appobject$N.reg
# Latent <- appobject$latent.trend
# rho.harv <- appobject$rho.harv
# rho.hunt <- appobject$rho.hunt
# n.species <- nrow(H.reg[[1]])
# 
# mu.N <- mu.H <- array(NA, dim = c(n.species, 2, n.year.f, n.samples))
# N.Est <- N.LCL <- N.UCL <- H.Est <- H.LCL <- H.UCL <- array(NA, dim = c(n.species, 2, n.year.f))
# 
# for(r in 1:2){
#     for(s in 1:n.species){
#       for(t in start.year.f:last.year.f){
#         mu.H[s,r,t-y.b,] <- rnorm(n.samples, H.reg[[r]]$A.hunt[s], H.reg[[r]]$A.hunt.sd[s]) +
#           econ.var*rnorm(n.samples, H.reg[[r]]$B.econ.hunt[s], H.reg[[r]]$B.econ.hunt.sd[s]) +
#           rnorm(n.samples, as.numeric(Latent$mean[[r]][s,t]), as.numeric(Latent$sd[[r]][s,t]))
#         
#         mu.N[s,r,t-y.b,] <- rnorm(n.samples, N.reg[[r]]$A.harv[s], N.reg[[r]]$A.harv.sd[s]) +
#           rnorm(n.samples, N.reg[[r]]$B.hunter[s], N.reg[[r]]$B.hunter.sd[s])*rnorm(n.samples, as.numeric(Latent$mean[[r]][s,t]), as.numeric(Latent$sd[[r]][s,t]))+
#           winter.var * rnorm(n.samples, N.reg[[r]]$B.ws.harv[s], N.reg[[r]]$B.ws.harv.sd[s]) +
#           drought.var * rnorm(n.samples, N.reg[[r]]$B.pdsi.harv[s], N.reg[[r]]$B.pdsi.harv.sd[s]) +
#           predator.var * rnorm(n.samples, N.reg[[r]]$B.bbs.harv[s], N.reg[[r]]$B.bbs.harv.sd[s])
#       }
#       for(t in 1:n.year.f){
#         H.Est[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(.5)))
#         H.LCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(.075)))
#         H.UCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(.935)))
#         
#         N.Est[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(.5)))
#         N.LCL[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(.075)))
#         N.UCL[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(.935)))
#       }
#     }
# }
# 
# BPH.Est <- N.Est/H.Est
# BPH.LCL <- N.LCL/H.LCL
# BPH.UCL <- N.UCL/H.UCL 
# 
# sem_out_comb <- function(f.arr){
#   as.data.frame(wrap(f.arr, map=list(NA, 2))) %>%
#     mutate(Species = rep(1:n.species, n.year.f),
#            Year = rep(start.year.f:last.year.f, each = n.species)) %>%
#     pivot_longer(cols = 1:2, names_to = "Region", names_prefix = "V") %>%
#     as.data.frame()
# }
# 
# N.df1 <- sem_out_comb(N.Est) %>% dplyr::rename(Est = value)
# N.df2 <- sem_out_comb(N.LCL) %>% dplyr::rename(LCL = value) %>%
#   merge(N.df1, ., by = c("Species", "Region", "Year"))
# N.df <- sem_out_comb(N.UCL) %>% dplyr::rename(UCL = value) %>%
#   merge(N.df2, ., by = c("Species", "Region", "Year"))
# 
# H.df1 <- sem_out_comb(H.Est) %>% dplyr::rename(Est = value)
# H.df2 <- sem_out_comb(H.LCL) %>% dplyr::rename(LCL = value) %>%
#   merge(H.df1, ., by = c("Species", "Region", "Year"))
# H.df <- sem_out_comb(H.UCL) %>% dplyr::rename(UCL = value) %>%
#   merge(H.df2, ., by = c("Species", "Region", "Year"))
# 
# BPH.df1 <- sem_out_comb(BPH.Est) %>% dplyr::rename(Est = value)
# BPH.df2 <- sem_out_comb(BPH.LCL) %>% dplyr::rename(LCL = value) %>%
#   merge(BPH.df1, ., by = c("Species", "Region", "Year"))
# BPH.df <- sem_out_comb(BPH.UCL) %>% dplyr::rename(UCL = value) %>%
#   merge(BPH.df2, ., by = c("Species", "Region", "Year"))
# 



# #Estimate
# mu.hunt.est <- Hunt.reg.coef %>%
#   merge(.,expand.grid(Species = 1:5, Region = 1:2, Year = 1:n.year.f),
#         by = c("Species", "Region"), all = T) %>%
#   mutate(Row = row_number()) %>%
#   rowwise() %>%
#   mutate(SPL = inprod(.[Row,3:14], Z.hunt[length(1976:(start.year.f-1))+Year,1:12,Species,Region])) %>%
#   mutate(Mu.hunt = A.hunt + B.econ.hunt*econ.var + SPL) %>%
#   as.data.frame() %>%
#   dplyr::select(Species, Region, Year, Mu.hunt) %>%
#   pivot_wider(names_from = Year, values_from = Mu.hunt, names_prefix = "Year", names_sep = "_")  %>%
#   group_split(Region, .keep = F)
# 
# 
# H.est <- list(as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)),
#               as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)))
# for(r in 1:2){
#   for(t in 1:n.year){
#     H.est[[r]][,t] <- as.numeric(unlist(exp(rmvnorm(1, mean = as.numeric(unlist(mu.hunt.est[[r]][,t+1])), sigma = as.matrix(rho.hunt.list[[r]])))))
#   }
# }
# 
# forecast.plot.input <- H.est %>% bind_rows() %>% mutate(Species = rep(1:5, 2), Region = rep(1:2, each = 5)) %>%
#   pivot_longer(cols = 1:n.year.f, names_prefix= "V") %>%
#   rename(Year = name, H.est = value) %>%
#   as.data.frame()
# 
# #Confidence Limits
# mu.hunt.est <- Hunt.reg.coef %>%
#   merge(.,expand.grid(Species = 1:5, Region = 1:2, Year = 1:n.year.f),
#         by = c("Species", "Region"), all = T) %>%
#   mutate(Row = row_number()) %>%
#   rowwise() %>%
#   mutate(SPL = inprod(.[Row,3:14], Z.hunt[length(1976:(start.year.f-1))+Year,1:12,Species,Region])) %>%
#   mutate(LCL.hunt = A.hunt + (B.econ.hunt*econ.var + SPL)) %>%
#   as.data.frame() %>%
#   dplyr::select(Species, Region, Year, Mu.hunt) %>%
#   pivot_wider(names_from = Year, values_from = Mu.hunt, names_prefix = "Year", names_sep = "_")  %>%
#   group_split(Region, .keep = F)
# 
# 
# H.est <- list(as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)),
#               as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)))
# for(r in 1:2){
#   for(t in 1:n.year){
#     H.est[[r]][,t] <- as.numeric(unlist(exp(rmvnorm(1, mean = as.numeric(unlist(mu.hunt.est[[r]][,t+1])), sigma = as.matrix(rho.hunt.list[[r]])))))
#   }
# }
# 
# H.est %>% bind_rows() %>% mutate(Species = rep(1:5, 2), Region = rep(1:2, each = 5)) %>%
#   pivot_longer(cols = 1:n.year.f, names_prefix= "V") %>%
#   rename(Year = name, H.est = value) %>%
#   as.data.frame()