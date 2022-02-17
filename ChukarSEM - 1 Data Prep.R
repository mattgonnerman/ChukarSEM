lapply(c("dplyr", "ggplot2", "reshape2", "reshape"), require, character.only = T) 

setwd('E:/Maine Drive/Projects/Gibson Chukkar/Code/PDSI')

temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)

CC    <- myfiles[[1]]
ELKO  <- myfiles[[2]]
ELY   <- myfiles[[3]]
EUR   <- myfiles[[4]]
FAL   <- myfiles[[5]]
RL    <- myfiles[[6]]
SS    <- myfiles[[7]]
TO    <- myfiles[[8]]
WIN   <- myfiles[[9]]

CC  $Station <- rep('CC')
ELKO$Station <- rep('ELKO')
ELY $Station <- rep('ELY')
EUR $Station <- rep('EUR')
FAL $Station <- rep('FAL')
RL  $Station <- rep('RL')
SS  $Station <- rep('SS')
TO  $Station <- rep('TO')
WIN $Station <- rep('WIN')



df <- rbind(CC,ELKO,ELY,EUR,FAL,RL,SS ,TO ,WIN )

df$TimeStep <- NULL

df$Year <- format(as.Date(df$Date, format="%m/%d/%Y"),"%Y")
df$Month <- format(as.Date(df$Date, format="%m/%d/%Y"),"%m")


df$Region <- recode(df$Station, 'CC' = 'West','ELKO' = 'East','ELY' = 'East','EUR' = 'East','FAL' = 'West','RL' = 'East','SS' = 'East','TO' = 'West','WIN' = 'West')
df$Season <- recode(df$Month,`01` = 'Winter',`02` = 'Winter',`03` = 'Winter',`04` = 'Early',`05` = 'Early',`06` = 'Early',`07` = 'Late',`08` = 'Late',`09` = 'Late',`10` = 'Winter',`11` = 'Winter',`12` = 'Winter' )

df$F.Year <- ifelse(df$Month == '10' | df$Month == '11' | df$Month == '12', as.numeric(as.character(df$Year)) + 1,as.numeric(as.character(df$Year)))


df_agg <- df %>% 
  group_by(Station,F.Year,Season) %>%
  summarize(m.PDSI = mean(PDSI, na.rm = TRUE))


df_cast <- cast(df_agg, F.Year + Season ~ Station, value = 'm.PDSI')

ggplot(df_agg, aes(x = F.Year, y = m.PDSI, color = Season)) +
  geom_point()+
  facet_wrap(~Station, ncol = 3) + theme_bw()

setwd('E:/Maine Drive/Projects/Gibson Chukkar/Code')

survey_data <- read.csv('Chukar_Surveys.csv')
harvest_data <- read.csv('all_species_harvest_data.csv')
harvest_data$Year <- as.numeric(as.character(harvest_data$Year))
harvest_data$Animals <- as.numeric(as.character(harvest_data$Animals))
harvest_data$Hunters <- as.numeric(as.character(harvest_data$Hunters))
harvest_data$H.Days <- as.numeric(as.character(harvest_data$H.Days  ))

lapply(c("jagsUI", "tidyverse", "nimble"), require, character.only = T) 

df <- gather(survey_data, "Population", "Count", 2:14)

site_list <- unique(df$Population)
year_list <- rep(1972:2017)

df_full <- data.frame("Population" = rep(site_list, times = length(year_list)), "Year" = rep(year_list, each = length(site_list)))

df_full <- left_join(df_full, df)

df_ch <- dcast(df_full, Population ~ Year)

ggplot() +
  geom_boxplot(data = df_full, aes(x = Year, y = Count, group = Year), fill = 'dodgerblue') +
  #geom_smooth(data = df_full, aes(x = Year, y = Count), fill = 'black',alpha = .5, color = 'dodgerblue1') +
  labs(y = 'Chukar covey survey counts') +
  theme_bw()

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

harvest_datae <- subset(harvest_data, grepl('Eastern', Section))

animal_wide <- cast(harvest_datae, Species_Code  ~ Year,  value = 'Animals', fun.aggregate = 'sum', na.rm = TRUE)
animal_wide[animal_wide==0] <- NA
species <- animal_wide$Species_Code 
animal_wide$Species_Code  <- NULL
animal_wide <- data.matrix(animal_wide)

hunters_wide <- cast(harvest_datae, Species_Code  ~ Year, value = 'Hunters', fun.aggregate = 'sum', na.rm = TRUE)
hunters_wide[hunters_wide==0] <-NA
hunters_wide$Species_Code  <- NULL
hunters_wide<- data.matrix(hunters_wide)

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

animal_widew <- animal_widew[-6,]
hunters_widew <- hunters_widew[-6, ]

require(abind)
upland <- abind(animal_widew,animal_wide, along = 3)
hunters <- abind(hunters_widew, hunters_wide, along = 3)

hunters <-hunters[-c(4,8),,]
upland <- upland[-c(4,8),,]

df <- cbind.data.frame(year = c(1976:2017), chukar = c(upland[2,,1]), hunter = c(hunters[2,,1])) 

ggplot(data = df) +
  geom_point(aes(x = year, y = log(chukar)), size = 3, shape = 21, fill = 'dodgerblue') +
  geom_point(aes(x = year, y = log(hunter)), size = 3, shape = 21, fill = 'goldenrod') +
  labs(y = 'Numbers of hunters and chukar (log-scaled)', x = 'Year') +
  theme_bw()

upland[7,10,] <- NA # Season closed
hunters[7,10,] <- NA # Season closed 

eastern_pdsi <- read.csv('Eastern_PDSIZ.csv')
eastern_pdsi <- subset(eastern_pdsi, Year > 1974)
winter_pdsie <- rep(NA, times = dim(eastern_pdsi)[1] - 1)

for(i in 1:length(winter_pdsie)){
  winter_pdsie[i] <- sum(eastern_pdsi[i, 12:13], eastern_pdsi[i+1,2:4] )/5
}

ebreeding_means <- data.frame(year = eastern_pdsi$Year, breeding_drought = rowMeans(eastern_pdsi[,5:9]))

eastern_pdsi <- subset(eastern_pdsi, Year > 1975)
ebreeding_means <- data.frame(year = eastern_pdsi$Year, breeding_drought = rowMeans(eastern_pdsi[,5:9]))

western_pdsi <- read.csv('Western_PDSIZ.csv')

western_pdsi <- subset(western_pdsi, Year > 1974)
winter_pdsiw <- rep(NA, times = dim(western_pdsi)[1] - 1)

for(i in 1:length(winter_pdsiw)){
  winter_pdsiw[i] <- sum(western_pdsi[i, 12:13], western_pdsi[i+1,2:4] )/5
}

western_pdsi <- subset(western_pdsi, Year > 1975)
wbreeding_means <- data.frame(year = western_pdsi$Year, breeding_drought = rowMeans(western_pdsi[,5:9]))

cut <- 42

pdsi  <- data.frame(east = ebreeding_means$breeding_drought[1:44],west = wbreeding_means$breeding_drought[1:44] )
wpdsi <- data.frame(east = winter_pdsie[1:(44+1)], west = winter_pdsiw[1:(44+1)])

#1976 - 2021
unemployment <- c(8.825,6.766666667,4.416666667,4.866666667,6.316666667,7.316666667,10.05,9.841666667,7.658333333,7.475,6.383333333,6.2,5.241666667,4.658333333,4.7,5.9,6.8,6.9,6.241666667,5.566666667,5.058333333,4.375,4.175,3.966666667,4.116666667,5.191666667,5.683333333,5.283333333,4.458333333,4.125,4.066666667,4.591666667,6.883333333,11.725,13.73333333,13.31666667,11.60833333,9.966666667,8.158333333,6.85,5.808333333,5.016666667,4.408333333,3.9,13.01666667,8.333333333)
unemployment_cut <- unemployment[1:(44+1)]
une <- scale(unemployment_cut)

general <- read.csv(file = 'overall_hunting_data.csv')
general_nv <- subset(general, ST == 'NV' & Year > 1975)
res <- (general_nv$Resident.Licenses)


econ_data <- read.csv('economic_data.csv')
econ_data_sub <- subset(econ_data, Year > 1975)
PDI <- (econ_data_sub$Per.Capita.Personal.Disposable.Income/10000)
GAS <- econ_data_sub$Gas.Price..July.Oct.

mean_ratio <- 2.581635
sd_ratio <- 0.8894599

ratio <- PDI/GAS