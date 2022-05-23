### Check Model Outputs against observed values

if(drop.rabbit == "Y"){
  check.species <- species[-c(4,7,8)]
}else{
  check.species <- species[-c(4,8)]
}


#Load all species harvest data
obs.HarvestData <- read.csv('./Data/all_species_harvest_data.csv') %>%
  mutate(Year = as.numeric(as.character(Year)),
         Animals = as.numeric(as.character(Animals)),
         Hunters = as.numeric(as.character(Hunters)),
         H.Days = as.numeric(as.character(H.Days))) %>%
  mutate(Species = recode(Species, 
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
         `GAMBEL'S_QUAIL` = 'GAQU')) %>%
  filter(Species %in% check.species) %>%
  dplyr::rename(Obs.N = Animals, Obs.H = Hunters, Region = Section) %>%
  mutate(Obs.BPH = Obs.N/Obs.H) %>%
  dplyr::select(-State, -H.Days) %>%
  filter(Region != "Southern")  %>%
  arrange(Region, Species, Year)
  
est.N <- test.N %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.N = Estimate, Est.N.LCL = LCL, Est.N.UCL = UCL) %>%
  filter(!is.na(Species))

est.H <- test.H %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.H = Estimate, Est.H.LCL = LCL, Est.H.UCL = UCL) %>%
  filter(!is.na(Species))

est.BPH <- test.bph %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.BPH = Estimate, Est.BPH.LCL = LCL, Est.BPH.UCL = UCL) %>%
  filter(!is.na(Species))  

est.check.df <- merge(est.N, est.H, by = c("Region", "Species", "Year"), all = T)
est.check.df <- merge(est.check.df, est.BPH, by = c("Region", "Species", "Year"), all = T)
est.check.df <- merge(est.check.df, obs.HarvestData, by = c("Region", "Species", "Year"), all = T) %>%
  arrange(Species, Year, Region)
is.num <- sapply(est.check.df, is.numeric)
est.check.df[is.num] <- lapply(est.check.df[is.num], round, 2)


ggplot(data = est.check.df, aes(x = Year)) +
  geom_line(aes(y = Obs.N, color = Region)) +
  geom_point(aes(y = Est.N, color = Region)) +
  geom_errorbar(aes(ymin = Est.N.LCL, ymax = Est.N.UCL, color = Region)) +
  facet_wrap(vars(Species), scales = "free_y") 
