y.start.fore <- 2015
y.stop.fore <- 2022

###MAKE SURE THE BELOW MATCHES THE ORIGINAL MODEL RUN
### Run Initial Data Management
lapply(c("parallel", "coda", "MCMCvis"), require, character.only = T)
cutoff.y <- 2016 #Last year from which data will be used
final.y <- 2017 #Last year to predict

drop.rabbit <- "N" #N to keep rabbit in harvest data correlation models

n.add.y <- final.y - cutoff.y
source("./ChukSEM - Data Prep - Predict.R")

#BBS Index
bbs.fore <- read.csv("./Data/bbs_indices.csv") %>%
  mutate(raven = scale(raven)[,1],
         rthawk = scale(rthawk)[,1],
         nharrier = scale(nharrier)[,1],
         pfalcon = scale(pfalcon)[,1]
  ) %>%
  filter(Year >= y.start.fore & Year <= y.stop.fore)

predict.fore1 <- rbind(bbs.fore %>% dplyr::select(year = Year, raven, nharrier),
                      data.frame(year = (1+max(bbs.fore$Year)):y.stop.fore,
                                 raven = NA, nharrier = NA))

#Drought Index
eastern.pdsi.fore <- read.csv('./Data/Eastern_PDSIZ.csv') %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer, Winter) %>%
  filter(Year %in% y.start.fore:y.stop.fore) %>%
  mutate(Region = "Eastern")

western.pdsi.fore <- read.csv('./Data/Western_PDSIZ.csv') %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer, Winter) %>%
  filter(Year %in% y.start.fore:y.stop.fore) %>%
  mutate(Region = "Western")

pdsi.fore <- rbind(eastern.pdsi.fore, western.pdsi.fore) %>%
  dplyr::rename(year = Year, pdsi_sum = Summer, pdsi_win = Winter, region = Region)

predict.fore2 <- merge(predict.fore1, pdsi.fore, by = "year", all.x = T) %>%
  arrange(year, region)

#Unemployment
une.mean <- attr(scale(unemployment$Rate), "scaled:center")
une.sd <- attr(scale(unemployment$Rate), "scaled:scale")
unemployment.fore <- read.csv("./Data/NEvadaUnemploymentUSBL.csv", colClasses = c("integer", "character", rep("numeric", 2), rep("integer", 3), "numeric")) %>%
  filter(Period == "Sep") %>%
  filter(Year %in% y.start.fore:y.stop.fore) %>%
  mutate(Rate = unemployment/labor.force) %>%
  mutate(Rate.Z = ifelse(is.na(Rate), NA, (Rate - une.mean)/une.sd)) %>%
  dplyr::select(year = Year, une = Rate.Z)
predict.fore3 <- merge(predict.fore2, unemployment.fore, by = "year", all.x = T)


#Number of Resident Licenses Sold
res.mean <- attr(scale(general_nv$Resident.Licenses), "scaled:center")
res.sd <- attr(scale(general_nv$Resident.Licenses), "scaled:scale")
general_nv.fore <- subset(general, ST == 'NV' & Year %in% y.start.fore:y.stop.fore) %>%
  mutate(Resident.Licenses.Z = ifelse(is.na(Resident.Licenses), NA, (Resident.Licenses - res.mean)/res.sd)) %>%
  dplyr::select(year = Year, res = Resident.Licenses.Z)
predict.fore4 <- merge(predict.fore3, general_nv.fore, by = "year", all.x = T)

#
gas.march.fore <- read.csv('./Data/NV_Gas.csv') %>%
  mutate(Date = as.Date(M.Y, format = '%m/%d/%Y')) %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date)) %>%
  filter(Month == 5) %>%
  arrange(Year) %>%
  filter(Year %in% y.start.fore:y.stop.fore) %>%
  dplyr::select(Year, Gas.May = Dollars.per.Gallon)

econ_data.fore <- read.csv('./Data/economic_data.csv') %>%
  dplyr::rename(PDI = Per.Capita.Personal.Disposable.Income, Gas.Unk = Gas.Price..July.Oct.) %>%
  merge(., gas.march.fore, by = "Year", all = T) %>%
  filter(Year %in% y.start.fore:y.stop.fore) %>%
  mutate(PDI = PDI/10000)

predict.fore$PDI <- econ_data$PDI
predict.fore$GAS <- econ_data$Gas.May
