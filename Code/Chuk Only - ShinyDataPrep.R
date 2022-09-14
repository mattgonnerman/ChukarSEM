lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr", "tidyr", "mvtnorm", "R.utils"), require, character.only = T)


#Load model outputs to get beta and sigma estimates
load(file = "./Output/ChukOnly - Final/NDOW_ChukOnly_SEM_output.rdata")
mcmcList1 <- files[[1]]
mcmcList2 <- files[[2]]
code <- files[[3]]
data <- files[[4]]

appobject.CH <- list()

appobject.CH$county_order <- county_order <- files[[5]]
appobject.CH$county_reg <- county_reg <- files[[6]]
appobject.CH$N.data <- data$n.harv
appobject.CH$H.data <- data$n.hunt
appobject.CH$data.all <- data

n.counties<- dim(data$n.hunt)[1]

### Observed Hunter and Harvest Data For Graphs
appobject.CH$H <- H.input <- data$n.hunt
appobject.CH$N <- N.input <- data$n.harv


### Hunter Effort Regression Coefficients
Hunt.reg.coef1 <- MCMCsummary(mcmcList2, 'alpha.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])"))) %>%
  mutate(ID = "b.0") %>%
  dplyr::select(County, A.hunt = mean, A.hunt.sd = sd) 
appobject.CH$H.reg <- Hunt.reg.coef <- MCMCsummary(mcmcList2, 'beta.econ.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.ECON") %>%
  dplyr::select(County, B.econ.hunt = mean, B.econ.hunt.sd = sd) %>%
  merge(Hunt.reg.coef1, ., by = c("County"), all =T)

appobject.CH$latent.trend <- latent.trend <- MCMCsummary(mcmcList2, 'latent.trend') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Latent = mean, Latent.sd = sd)


### Total Harvest Regression Coefficients
Harv.reg.coef1 <- MCMCsummary(mcmcList2, 'alpha.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])"))) %>%
  mutate(ID = "alpha.harv") %>%
  dplyr::select(County, A.harv = mean, A.harv.sd = sd) 
  
Harv.reg.coef2 <- MCMCsummary(mcmcList2, 'beta.wintsev.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(ID = "b.WINT") %>%
  dplyr::select(B.ws.harv = mean, B.ws.harv.sd = sd) %>%
  merge(., data.frame(County = 1:length(county_reg))) %>%
  merge(Harv.reg.coef1, ., by = c("County"), all =T)

Harv.reg.coef3 <- MCMCsummary(mcmcList2, 'beta.bbs.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(ID = "b.BBS") %>%
  dplyr::select(B.bbs.harv = mean, B.bbs.harv.sd = sd) %>%
  merge(., data.frame(County = 1:length(county_reg))) %>%
  merge(Harv.reg.coef2, ., by = c("County"), all =T)
  
appobject.CH$N.reg <- Harv.reg.coef <- MCMCsummary(mcmcList2, 'beta.hunter.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])"))) %>%
  dplyr::select(County, B.hunter = mean, B.hunter.sd = sd) %>%
  merge(Harv.reg.coef3, ., by = c("County"), all =T)

### SEM Predictor Values
appobject.CH$pred.econ <- MCMCsummary(mcmcList2, 'pred.econ.prime')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Year = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "econ.input") %>%
  dplyr::select(Year, Pred.Econ = mean, Pred.Econ.sd = sd)
appobject.CH$pred.bbs <- MCMCsummary(mcmcList2, 'pred.bbs.prime')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Year = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "econ.input") %>%
  dplyr::select(Year, Pred.Econ = mean, Pred.Econ.sd = sd)



save(appobject.CH, file = "./NDOWSEM/www/NDOW_ChukOnly_app_objects.R")
