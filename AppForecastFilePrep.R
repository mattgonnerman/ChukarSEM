lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr", "tidyr"), require, character.only = T)


#Load model outputs to get beta and sigma estimates
load(file = "./www/NDOW_Upland_SEM_output.rdata")
beta.output <- files[[2]]
rho.output <- files[[1]]
data.input <- files[[4]]

appobject <- list()

### Hunter Effort Covariance Matrix
appobject$rho.harv <- rho.harv.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 5, ncol = 5))
rho.harv.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.harv.est)){
  rho.harv.list[[rho.harv.est$Region[i]]][rho.harv.est$Species1[i], rho.harv.est$Species2[i]] <- rho.harv.est$Estimate[i]
}

### Hunter Effort Covariance Matrix
appobject$rho.hunt <- rho.hunt.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 5, ncol = 5))
rho.hunt.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.hunt.est)){
  rho.hunt.list[[rho.hunt.est$Region[i]]][rho.hunt.est$Species1[i], rho.hunt.est$Species2[i]] <- rho.hunt.est$Estimate[i]
}

### Observed Hunter and Harvest Data For Graphs
appobject$H <- H.input <- data.input$n.hunt
appobject$N <- N.input <- data.input$n.harv


### Spline Information
appobject$Z.harv <- Z.harv <- data.input$Z.harv
appobject$Z.hunt <- Z.hunt <- data.input$Z.hunt

### Hunter Effort Regression Coefficients
Hunt.reg.coef1 <- MCMCsummary(beta.output, 'beta.spl.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         K = sub('.*\\,', '', RowID)) %>%
  mutate(K = as.numeric(str_sub(K,1,nchar(K)-1)),
         ID = "b.SPL") %>%
  dplyr::select(Species, Region, K, mean) %>%
  pivot_wider(names_from = K, names_prefix = "Spl", values_from = c(mean)) 

Hunt.reg.coef2 <- MCMCsummary(beta.output, 'alpha.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID))%>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "b.0") %>%
  dplyr::select(Species, Region, A.hunt = mean, A.hunt.sd = sd) %>%
  merge(Hunt.reg.coef1, ., by = c("Species", "Region"), all =T)
  
appobject$H.reg <- Hunt.reg.coef <- MCMCsummary(beta.output, 'beta.econ.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.ECON") %>%
  dplyr::select(Species, B.econ.hunt = mean, B.econ.hunt.sd = sd) %>%
  merge(., expand.grid(Species = 1:5, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Hunt.reg.coef2, ., by = c("Species", "Region"), all =T)



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
  
appobject$N.reg <- Harv.reg.coef <- MCMCsummary(beta.output, 'beta.spl.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         K = sub('.*\\,', '', RowID)) %>%
  mutate(K = as.numeric(str_sub(K,1,nchar(K)-1))) %>%
  dplyr::select(Species, Region, K, mean) %>%
  pivot_wider(names_from = K, names_prefix = "Spl", values_from = c(mean)) %>%
  merge(Harv.reg.coef4, ., by = c("Species", "Region"), all =T)


save(appobject, file = "./www/NDOW_SEM_app_objects.R")

#### For Integration into Shiny
#These will represent the slider inputs. Will need to be reactive
require(nimble)
require(mvtnorm)

econ.var <- 0
winter.var <- 0
drought.var <- 0
predator.var <- 0
start.year.f <- 2018
last.year.f <- 2020
n.year.f <- length(start.year.f:last.year.f)

### Hunter Effort
#Estimate
mu.hunt.est <- Hunt.reg.coef %>%
  merge(.,expand.grid(Species = 1:5, Region = 1:2, Year = 1:n.year.f),
        by = c("Species", "Region"), all = T) %>%
  mutate(Row = row_number()) %>%
  rowwise() %>%
  mutate(SPL = inprod(.[Row,3:14], Z.hunt[length(1976:(start.year.f-1))+Year,1:12,Species,Region])) %>%
  mutate(Mu.hunt = A.hunt + B.econ.hunt*econ.var + SPL) %>%
  as.data.frame() %>%
  dplyr::select(Species, Region, Year, Mu.hunt) %>%
  pivot_wider(names_from = Year, values_from = Mu.hunt, names_prefix = "Year", names_sep = "_")  %>%
  group_split(Region, .keep = F)
  

H.est <- list(as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)),
              as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)))
for(r in 1:2){
  for(t in 1:n.year){
    H.est[[r]][,t] <- as.numeric(unlist(exp(rmvnorm(1, mean = as.numeric(unlist(mu.hunt.est[[r]][,t+1])), sigma = as.matrix(rho.hunt.list[[r]])))))
  }
}

forecast.plot.input <- H.est %>% bind_rows() %>% mutate(Species = rep(1:5, 2), Region = rep(1:2, each = 5)) %>%
  pivot_longer(cols = 1:n.year.f, names_prefix= "V") %>%
  rename(Year = name, H.est = value) %>%
  as.data.frame()

#Confidence Limits
mu.hunt.est <- Hunt.reg.coef %>%
  merge(.,expand.grid(Species = 1:5, Region = 1:2, Year = 1:n.year.f),
        by = c("Species", "Region"), all = T) %>%
  mutate(Row = row_number()) %>%
  rowwise() %>%
  mutate(SPL = inprod(.[Row,3:14], Z.hunt[length(1976:(start.year.f-1))+Year,1:12,Species,Region])) %>%
  mutate(LCL.hunt = A.hunt + (B.econ.hunt*econ.var + SPL) %>%
  as.data.frame() %>%
  dplyr::select(Species, Region, Year, Mu.hunt) %>%
  pivot_wider(names_from = Year, values_from = Mu.hunt, names_prefix = "Year", names_sep = "_")  %>%
  group_split(Region, .keep = F)


H.est <- list(as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)),
              as.data.frame(matrix(NA, nrow = 5, ncol = n.year.f)))
for(r in 1:2){
  for(t in 1:n.year){
    H.est[[r]][,t] <- as.numeric(unlist(exp(rmvnorm(1, mean = as.numeric(unlist(mu.hunt.est[[r]][,t+1])), sigma = as.matrix(rho.hunt.list[[r]])))))
  }
}

H.est %>% bind_rows() %>% mutate(Species = rep(1:5, 2), Region = rep(1:2, each = 5)) %>%
  pivot_longer(cols = 1:n.year.f, names_prefix= "V") %>%
  rename(Year = name, H.est = value) %>%
  as.data.frame()
