lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr", "tidyr", "mvtnorm", "R.utils"), require, character.only = T)


#Load model outputs to get beta and sigma estimates
load(file = "./Output/NDOW_Upland_SEM_output.rdata")
beta.output <- files[[2]]
rho.output <- files[[1]]
data.input <- files[[4]]

appobject <- list()

appobject$N.data <- data.input$n.harv
appobject$H.data <- data.input$n.hunt
appobject$data.all <- data.input

n.species<- dim(data.input$n.hunt)[1]

### Hunter Effort Covariance Matrix
rho.harv.est <- MCMCsummary(rho.output, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(rho.output, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = n.species, ncol = n.species))
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

blank.cor.df <- as.data.frame(matrix(NA, nrow = n.species, ncol = n.species))
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
  merge(., expand.grid(Species = 1:n.species, Region = 1:2), by = "Species", all.x = T) %>%
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
  merge(., expand.grid(Species = 1:n.species, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef1, ., by = c("Species", "Region"), all =T)

Harv.reg.coef3 <- MCMCsummary(beta.output, 'beta.pdsi.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.PDSI") %>%
  dplyr::select(Species, B.pdsi.harv = mean, B.pdsi.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:n.species, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef2, ., by = c("Species", "Region"), all =T)

Harv.reg.coef4 <- MCMCsummary(beta.output, 'beta.bbs.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.BBS") %>%
  dplyr::select(Species, B.bbs.harv = mean, B.bbs.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:n.species, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef3, ., by = c("Species", "Region"), all =T)
  
appobject$N.reg <- Harv.reg.coef <- MCMCsummary(beta.output, 'beta.hunter.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(Species, Region, B.hunter = mean, B.hunter.sd = sd) %>%
  merge(Harv.reg.coef4, ., by = c("Species", "Region"), all =T) %>%
  group_split(Region, .keep = F)


save(appobject, file = "./NDOWSEM/www/NDOW_SEM_app_objects.R")