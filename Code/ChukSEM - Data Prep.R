#####################################################################################
### All Species Harvest Data ###
#Load all species harvest data
harvest_data1 <- read.csv('./Data/all_species_harvest_data.csv') %>%
  mutate(Year = as.numeric(as.character(Year)),
         Animals = as.numeric(as.character(Animals)),
         Hunters = as.numeric(as.character(Hunters)),
         H.Days = as.numeric(as.character(H.Days))
  ) %>%
  filter(Year <= cutoff.y) %>%
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
  filter(Species %in% c("SAGR", "CHUK", "BLGR", "CAQU", "HUPA", "PHEA")) %>%
  dplyr::select(Section, Species, Year, N = Animals, H = Hunters)

list.all.y <- min(harvest_data1$Year):final.y
'%notin%' <- Negate('%in%')

harvest_data <- rbind(harvest_data1, expand.grid(Section = unique(harvest_data1$Section),
                                                 Species = unique(harvest_data1$Species),
                                                 Year = list.all.y[which(list.all.y %notin% unique(harvest_data1$Year))],
                                                 N = NA,
                                                 H = NA)) %>%
  arrange(Species, Section, Year)
  

#Subset data from Eastern section
harvest_datae <- subset(harvest_data, grepl('Eastern', Section))
#Reformat to number of individuals from each species harvested each year
animal_wide <- cast(harvest_datae, Species  ~ Year,  value = 'N', fun.aggregate = 'sum', na.rm = TRUE)
#Change 0 values to NA
animal_wide[animal_wide==0] <- NA
#Extract vector of species codes
species <- animal_wide$Species 
#Remove species code column from df
animal_wide$Species  <- NULL
#convert df to matrix
animal_wide <- data.matrix(animal_wide)

#Same thing as above but for number of hunters
hunters_wide <- cast(harvest_datae, Species  ~ Year, value = 'H', fun.aggregate = 'sum', na.rm = TRUE)
hunters_wide[hunters_wide==0] <-NA
hunters_wide$Species  <- NULL
hunters_wide<- data.matrix(hunters_wide)

#Same thing as previous but for western section
harvest_dataw <- subset(harvest_data, grepl('Western', Section))

animal_widew <- cast(harvest_dataw, Species  ~ Year,  value = 'N', fun.aggregate = 'sum', na.rm = TRUE)
animal_widew[animal_widew==0] <- NA
speciesw <- animal_widew$Species
animal_widew$Species  <- NULL
animal_widew <- data.matrix(animal_widew)

hunters_widew <- cast(harvest_dataw, Species  ~ Year, value = 'H', fun.aggregate = 'sum', na.rm = TRUE)
hunters_widew[hunters_widew==0] <-NA
hunters_widew$Species  <- NULL
hunters_widew<- data.matrix(hunters_widew)

#Combine datasets from East and West sections into arrays
upland <- abind(animal_wide, animal_widew,along = 3)
hunters <- abind(hunters_wide, hunters_widew, along = 3)

# #Add NA values to predict 
# upland <- abind(upland, array(NA, dim = c(nrow(upland), n.add.y, 2)), along = 2)
# hunters <- abind(hunters, array(NA, dim = c(nrow(hunters), n.add.y, 2)), along = 2)


#Birds per hunter from survey to inform missing years
surveybph <- read.csv("./Data/BirdsperHunterMissingYears.csv") %>%
  mutate(Species = ordered(Species, levels = species)) %>%
  arrange(Species) %>%
  filter(!is.na(Year)) %>%
  dplyr::select(Species, Year, BperH) %>%
  pivot_wider(names_from = "Year", values_from = "BperH") %>%
  dplyr::select(-Species) %>% as.matrix()

##################################################################################
### Chukar Survey Data ###
#Load chukar survey data
survey_data <- read.csv('./Data/Chukar_Surveys_data.csv')

#Number of surveys for each population per year
df <- gather(survey_data, "Population", "Count", 2:14)
# Where chukar surveys occurred
site_list <- unique(df$Population)
# Years to assess
year_list <- rep(1990:cutoff.y.chuk)
#Create new df for all population*year combinations
df_full <- data.frame("Population" = rep(site_list, times = length(year_list)), "Year" = rep(year_list, each = length(site_list)))
#Join associated number of surveys at each location in each year
df_full <- left_join(df_full, df)
#Format wider
df_ch <- dcast(df_full, Population ~ Year)
chukar <- df_ch %>% dplyr::select(as.character(year_list))
if(n.add.y >1){
  chukar <- cbind(chukar , as.data.frame(matrix(NA, nrow = nrow(chukar), ncol = n.add.y -1)))
}

chuk_site_data <- read.csv("./Data/Chukar_Surveys_locations.csv") %>%
  select(Population = Survey.Location, Region = NDOWREGION) %>%
  arrange(Population) %>%
  mutate(Code = ifelse(Region == "Western", 2, 1))
chuk.reg <- chuk_site_data[,3]


#####################################################################################
### Palmer Drought Severity Index (PDSI) ###
#Load year/month specific Eastern PDSI values (Z standardized?)
eastern_pdsi <- read.csv('./Data/Eastern_PDSIZ.csv') %>%
  filter(Year > 1975) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer) %>%
  mutate(Region = "Eastern")

western_pdsi <- read.csv('./Data/Western_PDSIZ.csv') %>%
  filter(Year > 1975) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer) %>%
  mutate(Region = "Western")

pdsi_df <- rbind(eastern_pdsi, western_pdsi) %>%
  pivot_wider(names_from = Region, values_from = c(Summer)) %>% 
  filter(Year <= final.y) %>%
  rbind(., expand.grid(Eastern = NA, Western = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  dplyr::select(-Year) %>%
  as.matrix(.)



#####################################################################################
### Economic Metrics ###
#Unemployment data, 1976 - 2021 (Bureau of Labor, https://www.bls.gov/regions/west/nevada.htm#eag)
unemployment <- read.csv("./Data/NEvadaUnemploymentUSBL.csv", colClasses = c("integer", "character", rep("numeric", 2), rep("integer", 3), "numeric")) %>%
  filter(Period == "Sep") %>%
  filter(Year > 1975) %>%
  mutate(Rate = unemployment/labor.force) %>%
  dplyr::select(Year, Une = Rate) %>%
  mutate(Une = scale(Une)[,1]) %>%
  rbind(., expand.grid(Une = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))]))

#Number of Resident Licenses Sold
general <- read.csv(file = './Data/overall_hunting_data.csv') %>%
  filter(ST == 'NV' & Year > 1975) %>%
  dplyr::select(Year, Licenses = Resident.Licenses) %>%
  mutate(Licenses = scale(Licenses)[,1]) %>%
  rbind(., expand.grid(Licenses = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))]))

#Gas Prices
#https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=EMA_EPM0_PWG_SNV_DPG&f=M
gas.march <- read.csv('./Data/NV_Gas.csv') %>%
  mutate(Date = as.Date(M.Y, format = '%m/%d/%Y')) %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date)) %>%
  filter(Month == 5) %>%
  arrange(Year) %>%
  dplyr::select(Year, Gas.May = Dollars.per.Gallon) %>%
  mutate(Gas.May = scale(Gas.May)[,1]) %>%
  rbind(., expand.grid(Gas.May = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))]))

#Personal Disposalable Income/Merge All
econ_data <- read.csv('./Data/economic_data.csv') %>%
  dplyr::select(Year, PDI = Per.Capita.Personal.Disposable.Income) %>%
  filter(Year >1975) %>%
  mutate(PDI = scale(PDI)[,1]) %>%
  rbind(., expand.grid(PDI = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  merge(., gas.march, by = "Year", all = T)%>%
  merge(., general, by = "Year", all = T) %>%
  merge(., unemployment, by = "Year", all = T) %>%
  filter(Year <= final.y) %>%
  dplyr::select(Gas.May, Une, Licenses, PDI)


#####################################################################################
#Winter Severity (AWSSI)
awssi.df <- read.csv("./Data/Nevada AWSSI.csv") %>%
  select(Station, Region, Start, End, AWSSI) %>%
  mutate(Start = as.Date(Start, "%m/%d/%Y")) %>%
  mutate(Year = lubridate::year(Start)) %>%
  group_by(Region, Year) %>%
  filter(!is.na(AWSSI)) %>%
  summarise(AWSSI = mean(AWSSI)) %>%
  filter(Year > 1974 & Region != "Southern") %>%
  mutate(AWSSI = scale(AWSSI)[,1]) %>%
  pivot_wider(names_from = Region, values_from = AWSSI) %>%
  rbind(., expand.grid(Eastern = NA, Western = NA, 
                       Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  filter(Year <= final.y -1) %>%
  dplyr::select(-Year) %>%
  as.matrix(.)


#####################################################################################
#BBS Data
bbs.df <- read.csv("./Data/bbs_indices.csv") %>%
  mutate(raven = scale(raven)[,1],
         rthawk = scale(rthawk)[,1],
         nharrier = scale(nharrier)[,1],
         pfalcon = scale(pfalcon)[,1]
         ) %>%
  rbind(., expand.grid(raven = NA, rthawk = NA, nharrier = NA, pfalcon = NA, 
                       Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  filter(Year <= final.y -1) %>%
  dplyr::select(-Year) %>%
  as.matrix()


#####################################################################################
#Bobcat Productivity Data
bobcat.df <- read.csv("./Data/NV Bobcat Prod.csv") %>%
  dplyr::rename(bobcat = Kper100F) %>% dplyr::select(-Period) %>%
  mutate(bobcat = scale(bobcat)[,1]) %>%
  rbind(., expand.grid(bobcat = NA, Year = c(min(list.all.y) -1, list.all.y)[c(1, 1+which(list.all.y %notin% unique(.$Year)))])) %>%
  filter(Year <= final.y -1) %>%
  arrange(Year)

