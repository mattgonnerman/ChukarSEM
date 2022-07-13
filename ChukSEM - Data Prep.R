# Load necessary packages
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "LaplacesDemon"),
       require, character.only = T) 

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
### All Species Harvest Data ###
#Load all species harvest data
harvest_data <- read.csv('./Data/all_species_harvest_data.csv') %>%
  mutate(Year = as.numeric(as.character(Year)),
         Animals = as.numeric(as.character(Animals)),
         Hunters = as.numeric(as.character(Hunters)),
         H.Days = as.numeric(as.character(H.Days))
  ) %>%
  filter(Year <= cutoff.y)

#Add species codes
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

#Subset data from Eastern section
harvest_datae <- subset(harvest_data, grepl('Eastern', Section))
#Reformat to number of individuals from each species harvested each year
animal_wide <- cast(harvest_datae, Species_Code  ~ Year,  value = 'Animals', fun.aggregate = 'sum', na.rm = TRUE)
#Change 0 values to NA
animal_wide[animal_wide==0] <- NA
#Extract vector of species codes
species <- animal_wide$Species_Code 
#Remove species code column from df
animal_wide$Species_Code  <- NULL
#convert df to matrix
animal_wide <- data.matrix(animal_wide)

#Same thing as above but for number of hunters
hunters_wide <- cast(harvest_datae, Species_Code  ~ Year, value = 'Hunters', fun.aggregate = 'sum', na.rm = TRUE)
hunters_wide[hunters_wide==0] <-NA
hunters_wide$Species_Code  <- NULL
hunters_wide<- data.matrix(hunters_wide)

#Same thing as previous but for western section
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

#Remove MOQU
animal_widew <- animal_widew[-6,]
hunters_widew <- hunters_widew[-6, ]

#Combine datasets from East and West sections into arrays
upland <- abind(animal_wide, animal_widew,along = 3)
hunters <- abind(hunters_wide, hunters_widew, along = 3)

#Remove additional rows associated with species we don't want to estimate
if(drop.rabbit == "Y"){
  hunters <-hunters[-c(4,7,8),,]
  upland <- upland[-c(4,7,8),,]
}else{
  hunters <-hunters[-c(4,8),,]
  upland <- upland[-c(4,8),,]
}


#Change closed season values to NA
upland[6,10,] <- NA # Season closed
hunters[6,10,] <- NA # Season closed 

#Add NA values to predict 
upland <- abind(upland, array(NA, dim = c(nrow(upland), n.add.y, 2)), along = 2)
hunters <- abind(hunters, array(NA, dim = c(nrow(hunters), n.add.y, 2)), along = 2)


#####################################################################################
### Sage Grouse Wing Bee
#Format WingBee Data
sg.wingb <- read.csv("./Data/SG_WingData_2004-2020.csv") %>%
  select(Region = NDOWREGION, Year, AHY.Male, AHY.Female, HY.Male, HY.Female) %>%
  group_by(Region, Year) %>%
  summarize(AHY.F = sum(AHY.Female, na.rm = T),
            HY.F = sum(HY.Female, na.rm = T),
            HY.M = sum(HY.Male, na.rm = T)) %>%
  mutate(HY = HY.M + HY.F) %>%
  select(Region, Year, HY, AHY.F) %>%
  pivot_wider(names_from = "Region", values_from = c("AHY.F", "HY")) %>%
  filter(Year <= cutoff.y)

wing.b.ahy <- t(sg.wingb[,2:3])
wing.b.ahy <- cbind(wing.b.ahy , matrix(NA, nrow = 2, ncol = n.add.y))
wing.b.hy <- t(sg.wingb[,4:5])
wing.b.hy <- cbind(wing.b.hy , matrix(NA, nrow = 2, ncol = n.add.y))
n.years.sg <- ncol(wing.b.ahy)
time.shift.sg <- min(sg.wingb$Year) - min(harvest_data$Year) - 1


#####################################################################################
### Palmer Drought Severity Index (PDSI) ###
#Load year/month specific Eastern PDSI values (Z standardized?)
eastern_pdsi <- read.csv('./Data/Eastern_PDSIZ.csv') %>%
  filter(Year > 1974) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer, Winter) %>%
  filter(Year %in% 1975:cutoff.y) %>%
  mutate(Region = "Eastern")

western_pdsi <- read.csv('./Data/Western_PDSIZ.csv') %>%
  filter(Year > 1974) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer, Winter) %>%
  filter(Year %in% 1975:cutoff.y) %>%
  mutate(Region = "Western")

pdsi_df <- rbind(eastern_pdsi, western_pdsi) %>%
  pivot_wider(names_from = Region, values_from = c(Winter, Summer))

#Average breeding and winter DPSI values for each section
pdsi <- as.matrix(pdsi_df %>% select(Summer_Eastern, Summer_Western)) 
pdsi <- rbind(pdsi, matrix(NA, ncol = 2, nrow = n.add.y))
wpdsi <- as.matrix(pdsi_df[-1,] %>% select(Winter_Eastern, Winter_Western))
wpdsi <- rbind(wpdsi, matrix(NA, ncol = 2, nrow = n.add.y))


#####################################################################################
### Economic Metrics ###
#Unemployment data, 1976 - 2021 (Bureau of Labor, https://www.bls.gov/regions/west/nevada.htm#eag)
unemployment <- read.csv("./Data/NEvadaUnemploymentUSBL.csv", colClasses = c("integer", "character", rep("numeric", 2), rep("integer", 3), "numeric")) %>%
  filter(Period == "Sep") %>%
  filter(Year > 1975) %>%
  filter(Year <= cutoff.y) %>%
  mutate(Rate = unemployment/labor.force)
une <- c(scale(unemployment$Rate)[,1], rep(NA, n.add.y))
# ggplot(data = data.frame(UNE= une, Year = 1:length(une)), aes(x = Year, y = UNE)) +
#   geom_line()

#Number of Resident Licenses Sold
general <- read.csv(file = './Data/overall_hunting_data.csv')
general_nv <- subset(general, ST == 'NV' & Year > 1975 & Year <= cutoff.y)
res <- c(scale(general_nv$Resident.Licenses)[,1], rep(NA, n.add.y))

#Gas Prices
#https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=EMA_EPM0_PWG_SNV_DPG&f=M
gas.march <- read.csv('./Data/NV_Gas.csv') %>%
  mutate(Date = as.Date(M.Y, format = '%m/%d/%Y')) %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date)) %>%
  filter(Month == 5) %>%
  arrange(Year) %>%
  filter(Year < (final.y + 1)) %>%
  dplyr::select(Year, Gas.May = Dollars.per.Gallon)

econ_data <- read.csv('./Data/economic_data.csv') %>%
  dplyr::rename(PDI = Per.Capita.Personal.Disposable.Income, Gas.Unk = Gas.Price..July.Oct.) %>%
  merge(., gas.march, by = "Year", all = T) %>%
  filter(Year >1975 & Year < (final.y + 1)) %>%
  mutate(PDI = PDI/10000)

PDI <- econ_data$PDI
GAS <- econ_data$Gas.May


#Objects not used but values are represented in model
mean_ratio <- 2.581635
sd_ratio <- 0.8894599

#Calculate ratio of disposable income to gas price
ratio <- PDI/GAS


#####################################################################################
#Winter Severity (AWSSI)
awssi.df <- read.csv("./Data/Nevada AWSSI.csv") %>%
  select(Station, Region, Start, End, AWSSI) %>%
  mutate(Start = as.Date(Start, "%m/%d/%Y")) %>%
  mutate(Year = lubridate::year(Start)) %>%
  group_by(Region, Year) %>%
  filter(!is.na(AWSSI)) %>%
  summarise(AWSSI = mean(AWSSI)) %>%
  filter(Year > 1974 & Year <= cutoff.y & Region != "Southern") %>%
  mutate(AWSSI = scale(AWSSI)[,1]) %>%
  pivot_wider(names_from = Region, values_from = AWSSI) 

awssi <- cbind(t(awssi.df %>%
             select(-Year)), matrix(NA, nrow = 2, ncol = n.add.y))
# awssi <- awssi[,-ncol(awssi)] #Only need to include if you don't want concurrent winter severity


#####################################################################################
#BBS Data
bbs.df <- read.csv("./Data/bbs_indices.csv") %>%
  mutate(raven = scale(raven)[,1],
         rthawk = scale(rthawk)[,1],
         nharrier = scale(nharrier)[,1],
         pfalcon = scale(pfalcon)[,1]
         ) %>%
  filter(Year <= cutoff.y) %>%
  as.matrix()

bbs.df <- as.data.frame(rbind(bbs.df, matrix(NA, ncol = 5, nrow = n.add.y)))
