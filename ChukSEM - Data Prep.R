# Load necessary packages
lapply(c("dplyr", "ggplot2", "reshape2", "reshape", "jagsUI", "tidyverse", "nimble",
         "abind", "LaplacesDemon"),
       require, character.only = T) 


### Chukar Survey Data ###
#Load chukar survey data
survey_data <- read.csv('./Data/Chukar_Surveys_data.csv')

#Number of surveys for each population per year
df <- gather(survey_data, "Population", "Count", 2:14)
# Where chukar surveys occurred
site_list <- unique(df$Population)
# Years to assess
year_list <- rep(1972:cutoff.y)
#Create new df for all population*year combinations
df_full <- data.frame("Population" = rep(site_list, times = length(year_list)), "Year" = rep(year_list, each = length(site_list)))
#Join associated number of surveys at each location in each year
df_full <- left_join(df_full, df)
#Format wider
df_ch <- dcast(df_full, Population ~ Year)
chukar <- df_ch[,-c(1:19)]
chukar <- cbind(chukar , matrix(NA, nrow = nrow(chukar), ncol = n.add.y))

chuk_site_data <- read.csv("./Data/Chukar_Surveys_locations.csv") %>%
  select(Population = Survey.Location, Region = NDOWREGION) %>%
  arrange(Population) %>%
  mutate(Code = ifelse(Region == "Western", 2, 1))
chuk.reg <- chuk_site_data[,3]

#Show distribution of surveys across years
# ggplot() +
#   geom_boxplot(data = df_full, aes(x = Year, y = Count, group = Year), fill = 'dodgerblue') +
#   #geom_smooth(data = df_full, aes(x = Year, y = Count), fill = 'black',alpha = .5, color = 'dodgerblue1') +
#   labs(y = 'Chukar covey survey counts') +
#   theme_bw()



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


### Palmer Drought Severity Index (PDSI) ###
#Load and format original data files
temp <- list.files(path = "./Data/PDSI/",
                   pattern="*.csv")
temp <- paste("./Data/PDSI/", temp, sep = "")
myfiles <- lapply(temp, read.csv)

CC    <- myfiles[[1]]
ELKO  <- myfiles[[2]]
ELY   <- myfiles[[3]]
EUR   <- myfiles[[4]]
FAL   <- myfiles[[5]]
RL    <- myfiles[[6]]
SS    <- myfiles[[7]]
TO    <- myfiles[[8]]
WIN   <- myfiles[[9]]

#Add station ID
CC  $Station <- rep('CC')
ELKO$Station <- rep('ELKO')
ELY $Station <- rep('ELY')
EUR $Station <- rep('EUR')
FAL $Station <- rep('FAL')
RL  $Station <- rep('RL')
SS  $Station <- rep('SS')
TO  $Station <- rep('TO')
WIN $Station <- rep('WIN')

#Combine into individual DPSI data into single DF
df <- rbind(CC,ELKO,ELY,EUR,FAL,RL,SS ,TO ,WIN )
#remove timestep column
df$TimeStep <- NULL 
#add year column
df$Year <- format(as.Date(df$Date, format="%m/%d/%Y"),"%Y")
#add month column
df$Month <- format(as.Date(df$Date, format="%m/%d/%Y"),"%m")
#Recode station ID to region information
df$Region <- recode(df$Station, 'CC' = 'West','ELKO' = 'East','ELY' = 'East','EUR' = 'East','FAL' = 'West','RL' = 'East','SS' = 'East','TO' = 'West','WIN' = 'West')
#Recode month to season and fiscal year information
df$Season <- recode(df$Month,`01` = 'Winter',`02` = 'Winter',`03` = 'Winter',`04` = 'Early',`05` = 'Early',`06` = 'Early',`07` = 'Late',`08` = 'Late',`09` = 'Late',`10` = 'Winter',`11` = 'Winter',`12` = 'Winter' )
df$F.Year <- ifelse(df$Month == '10' | df$Month == '11' | df$Month == '12', as.numeric(as.character(df$Year)) + 1,as.numeric(as.character(df$Year)))

#Simplify df to mean PDSI for each station/fiscal year/season combination
df_agg <- df %>% 
  group_by(Station,F.Year,Season) %>%
  summarize(m.PDSI = mean(PDSI, na.rm = TRUE))

#Reformat to wider table
df_cast <- cast(df_agg, F.Year + Season ~ Station, value = 'm.PDSI')

#Show distribution of data across years
# ggplot(df_agg, aes(x = F.Year, y = m.PDSI, color = Season)) +
#   geom_point()+
#   facet_wrap(~Station, ncol = 3) + theme_bw()

### Palmer Drought Severity Index (PDSI)
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

# #Load year/month specific Eastern PDSI values (Z standardized?)
# eastern_pdsi <- read.csv('./Data/Eastern_PDSIZ.csv')
# # Create vector of year-specific average winter PDSI values
# eastern_pdsi <- subset(eastern_pdsi, Year > 1974)
# winter_pdsie <- rep(NA, times = dim(eastern_pdsi)[1] - 1)
# for(i in 1:length(winter_pdsie)){
#   winter_pdsie[i] <- sum(eastern_pdsi[i, 12:13], eastern_pdsi[i+1,2:4] )/5
# }
# #Subset to years after 1975
# eastern_pdsi <- subset(eastern_pdsi, Year > 1975)
# #Create dataframe of year and mean DPSI during breeding season
# ebreeding_means <- data.frame(year = eastern_pdsi$Year, breeding_drought = rowMeans(eastern_pdsi[,5:9]))
# 
# #Same as above but for western section
# western_pdsi <- read.csv('./Data/Western_PDSIZ.csv')
# western_pdsi <- subset(western_pdsi, Year > 1974)
# winter_pdsiw <- rep(NA, times = dim(western_pdsi)[1] - 1)
# for(i in 1:length(winter_pdsiw)){
#   winter_pdsiw[i] <- sum(western_pdsi[i, 12:13], western_pdsi[i+1,2:4] )/5
# }
# western_pdsi <- subset(western_pdsi, Year > 1975)
# wbreeding_means <- data.frame(year = western_pdsi$Year, breeding_drought = rowMeans(western_pdsi[,5:9]))
# 
cut <- length(1976:cutoff.y) + n.add.y #Reference used to subset dataframes later
# 
# #Average breeding and winter DPSI values for each section
# pdsi  <- data.frame(east = ebreeding_means$breeding_drought[1:44],west = wbreeding_means$breeding_drought[1:44] )
# wpdsi <- data.frame(east = winter_pdsie[1:(44+1)], west = winter_pdsiw[1:(44+1)])

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

#Rabbit Harvest
rabbits <- as.matrix(read.csv('./Data/all_species_harvest_data.csv') %>%
                       filter(Species == "RABBIT", Section !="Southern") %>%
                       select(Region = Section, Year, Rabbits = Animals) %>%
                       filter(Year <= cutoff.y) %>%
                       mutate(Rabbits = scale(Rabbits)) %>%
                       pivot_wider(values_from = "Rabbits", names_from = "Region") %>%
                       select(Eastern, Western)) #row 1 = 1976
rabbits <- rbind(rabbits, matrix(NA, ncol = 2, nrow = n.add.y))

# ggplot(data = as.data.frame(rabbits) %>% mutate(Year = 1976:final.y), aes(x = Year)) +
#   geom_line(aes(y = Eastern), color = "red") +
#   geom_line(aes(y = Western), color = "blue")

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


# Hunter days (a covariate for total harvest model)
# harvest_data <- read.csv('./Data/all_species_harvest_data.csv') %>%
#   mutate(Year = as.numeric(as.character(Year)),
#          Animals = as.numeric(as.character(Animals)),
#          Hunters = as.numeric(as.character(Hunters)),
#          H.Days = as.numeric(as.character(H.Days))) %>%
#   filter(Year <= cutoff.y) %>% 
#   filter(Section == "Western") %>%
#   mutate(Species = recode(Species, 
#                           `SAGE_GROUSE` = "SAGR",
#                           `SOOTY_GROUSE` = "BLGR",
#                           `RUFFED_GROUSE` = "RUGR",
#                           `CHUKAR` = "CHUK",
#                           `HUNGARIAN_PART` = 'HUPA',
#                           `CALIFORNIA_QUAIL` = 'CAQU',
#                           `MOUNTAIN_QUAIL` = 'MOQU',
#                           `PHEASANTS` = 'PHEA',
#                           `RABBIT` = 'RABB',
#                           `DOVE` = 'DOVE',
#                           `DUSKY_GROUSE` = 'BLGR',
#                           `GAMBEL'S_QUAIL` = 'GAQU')) %>%
#   filter(Species %in% if(drop.rabbit == "Y"){species[-c(4,7,8)]}else{species[-c(4,8)]}) %>%
#   pivot_wider(id_cols = Species, names_from = Year, values_from = H.Days)
