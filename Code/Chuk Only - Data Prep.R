#####################################################################################
## Harvest Data
chukharv_data <- read.csv('./Data/ChukarHarvestData.csv') %>%
  mutate(N = round(N),
         H = round(H)) %>%
  filter(Year <= cutoff.y) %>%
  arrange(County)

county_order1 <- unique(chukharv_data$County)
# Regions, 2 == West, 1 == East, 3 == South
county_reg1 <- c(2, 2, 3, 2, 1, 3, 1, 2, 1, 3, 2, 2, 3, 2, 2, 2, 1)

chukharv_data <- rbind(chukharv_data,
                       expand.grid(County = county_order1, N = NA, H = NA,
                                   Year = list.all.y[which(list.all.y %notin% unique(chukharv_data$Year))])) %>%
  rowwise() %>%
  #Identify southern region and remove for now.
  mutate(Region = county_reg1[which(County == county_order1)]) %>%
  filter(Region != 3) %>% dplyr::select(-Region)

county_order <- unique(chukharv_data$County)
# Regions, 2 == West, 1 == East, 3 == South
county_reg <- c(2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 1)
n_county <- length(county_order)

chuk_hunt <- chukharv_data %>%
  arrange(County, Year) %>%
  dplyr::select(-N) %>%
  pivot_wider(values_from = H, names_from = Year) %>%
  mutate(C_Order = which(county_order == County)) %>%
  arrange(C_Order) %>%
  dplyr::select(-County, -C_Order) %>%
  as.data.frame() %>%
  mutate_all(as.numeric)

chuk_harv <- chukharv_data %>%
  arrange(County, Year) %>%
  dplyr::select(-H) %>%
  pivot_wider(values_from = N, names_from = Year) %>%
  mutate(C_Order = which(county_order == County)) %>%
  arrange(C_Order) %>%
  dplyr::select(-County, -C_Order) %>%
  as.data.frame() %>%
  mutate_all(as.numeric)

n.year.harv <- ncol(chuk_hunt)


#Birds per hunter from survey to inform missing years
surveybph.df <- read.csv("./Data/BirdsperHunterMissingYears.csv") %>%
  filter(!is.na(Year), Species == "CHUK") %>%
  dplyr::select(Species, Year, BperH) %>%
  pivot_wider(names_from = "Year", values_from = "BperH") %>%
  dplyr::select(-Species) %>% as.matrix()
surveybph <- as.numeric(surveybph.df[1,])

##################################################################################
### Chukar Survey Data ###
chuk_site_ID <- read.csv("./Data/Chukar_Surveys_locations.csv") %>%
  dplyr::select(Population = Survey.Location, County = COUNTY) %>%
  mutate(CountyID = match(County, county_order1, nomatch = NA)) %>%
  mutate(RegionID = county_reg1[CountyID])

chuksiteabun_data <- read.csv('./Data/Chukar_Surveys_data.csv') %>%
  gather("Population", "Count", 2:14) %>%
  filter(Year %in% 1990:cutoff.y.chuk) %>%
  rbind(., expand.grid(Population = unique(.$Population),
                       Year = list.all.y[which(list.all.y %notin% unique(.$Year))],
                       Count = NA)) %>%
  merge(., chuk_site_ID, by = "Population", all = T)  %>%
  filter(RegionID != 3)%>%
  dplyr::select(-RegionID, -CountyID, -County) %>%
  merge(.,expand.grid(Year = c(min(.$Year):max(.$Year)),
                      Population = unique(.$Population))) %>%
  pivot_wider(values_from = Count, names_from = Year) %>%
  merge(., chuk_site_ID, by = "Population", all.x = T) %>%
  mutate(CountyID = match(County, county_order, nomatch = NA)) %>%
  mutate(RegionID = county_reg[CountyID])

survey.abun <- as.matrix(chuksiteabun_data %>% dplyr::select(-Population, -County, -CountyID, -RegionID))
survey.county <- chuksiteabun_data$CountyID
survey.reg <- chuksiteabun_data$RegionID
n.year.surv <- ncol(survey.abun)


#####################################################################################
### Palmer Drought Severity Index (PDSI) ###
first.y.cov <- min(unique(chukharv_data$Year))

#Load year/month specific Eastern PDSI values (Z standardized?)
eastern_pdsi <- read.csv('./Data/Eastern_PDSIZ.csv') %>%
  filter(Year >= first.y.cov) %>%
  mutate(Summer = (Apr + May + Jun + Jul + Aug)/5,
         Winter = (Nov + Dec + lead(Jan) + lead(Feb) + lead(Mar))/5) %>%
  select(Year, Summer) %>%
  mutate(Region = "Eastern")

western_pdsi <- read.csv('./Data/Western_PDSIZ.csv') %>%
  filter(Year >= first.y.cov) %>%
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
  filter(Year >= first.y.cov) %>%
  mutate(Rate = unemployment/labor.force) %>%
  dplyr::select(Year, Une = Rate) %>%
  mutate(Une = scale(Une)[,1]) %>%
  rbind(., expand.grid(Une = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))]))

#Number of Resident Licenses Sold
general <- read.csv(file = './Data/overall_hunting_data.csv') %>%
  filter(ST == 'NV' & Year >= first.y.cov) %>%
  dplyr::select(Year, Licenses = Resident.Licenses) %>%
  mutate(Licenses = scale(Licenses)[,1]) %>%
  add_row(Year = (max(.$Year)+1):final.y, Licenses = rep(NA, length((max(.$Year)+1):final.y)))

#Gas Prices
#https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=EMA_EPM0_PWG_SNV_DPG&f=M
gas.march <- read.csv('./Data/NV_Gas.csv') %>%
  mutate(Date = as.Date(M.Y, format = '%m/%d/%Y')) %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date)) %>%
  filter(Month == 5) %>%
  filter(Year >= first.y.cov) %>%
  arrange(Year) %>%
  dplyr::select(Year, Gas.May = Dollars.per.Gallon) %>%
  mutate(Gas.May = scale(Gas.May)[,1]) %>%
  rbind(., expand.grid(Gas.May = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))]))

#Personal Disposalable Income/Merge All
econ_data <- read.csv('./Data/economic_data.csv') %>%
  dplyr::rename(PDI = Per.Capita.Personal.Disposable.Income) %>%
  filter(Year >= first.y.cov) %>%
  mutate(PDI = scale(PDI)[,1]) %>%
  dplyr::select(-3) %>%
  rbind(., expand.grid(PDI = NA, Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  merge(., gas.march, by = "Year", all = T)%>%
  merge(., general, by = "Year", all = T) %>%
  merge(., unemployment, by = "Year", all = T) %>%
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
  filter(Year >= first.y.cov & Region != "Southern") %>%
  mutate(AWSSI = scale(AWSSI)[,1]) %>%
  pivot_wider(names_from = Region, values_from = AWSSI) %>%
  rbind(., expand.grid(Eastern = NA, Western = NA, 
                       Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  dplyr::select(-Year) %>%
  as.matrix(.)


#####################################################################################
#BBS Data
bbs.df <- read.csv("./Data/bbs_indices.csv") %>%
  filter(Year >= first.y.cov) %>%
  mutate(raven = scale(raven)[,1],
         rthawk = scale(rthawk)[,1],
         nharrier = scale(nharrier)[,1],
         pfalcon = scale(pfalcon)[,1]) %>%
  rbind(., expand.grid(raven = NA, rthawk = NA, nharrier = NA, pfalcon = NA, 
                       Year = list.all.y[which(list.all.y %notin% unique(.$Year))])) %>%
  dplyr::select(-Year) %>%
  as.matrix()


#####################################################################################
#Bobcat Productivity Data
bobcat.df <- read.csv("./Data/NV Bobcat Prod.csv") %>%
  dplyr::rename(bobcat = Kper100F) %>% dplyr::select(-Period) %>%
  filter(Year >= min(list.all.y) -1) %>%
  mutate(bobcat = scale(bobcat)[,1]) %>%
  rbind(., expand.grid(bobcat = NA, Year = c(min(list.all.y) -1, list.all.y)[which(c(min(list.all.y) -1, list.all.y) %notin% unique(.$Year))])) %>%
  filter(Year <= final.y -1) %>%
  arrange(Year)


#####################################################################################
### Specify Additional Data Inputs
require(splines)
bs_bbase <- function(x, xl = min(x, na.rm = TRUE), xr = max(x, na.rm=TRUE), nseg = 10, deg = 3) {
  # Compute the length of the partitions
  dx <- (xr - xl) / nseg
  # Create equally spaced knots
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  # Use bs() function to generate the B-spline basis
  get_bs_matrix <- matrix(bs(x, knots = knots, degree = deg, Boundary.knots = c(knots[1], knots[length(knots)])), nrow = length(x))
  # Remove columns that contain zero only
  bs_matrix <- get_bs_matrix[, -c(1:deg, ncol(get_bs_matrix):(ncol(get_bs_matrix) - deg))]
  
  return(bs_matrix)
}

nseg <- 15 #Number of spline segments
time <- 1:cut
B <- bs_bbase(x = time, nseg = nseg)


#####################################################################################
### Specify Additional Constants Inputs
n.counties <- dim(chuk_hunt)[1]
sig = rgamma(n.counties,1,1)
Lambda <- Lambda2 <- diag(sig)
I2 <- I <- diag(n.counties) #identity matrix


#####################################################################################
### Specify Initial Values
## Hunter Effort
nu = n.counties + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.counties, n.counties)
for (j in 1:n.counties){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)

##Total Harvest
nu = n.counties + 1 #degrees of freedom for rinvwishart
Q = rinvwishart(nu, I)    # note output is an array
Delta = matrix(0, n.counties, n.counties)
for (j in 1:n.counties){
  Delta[j,j] = Q[j,j]^(-0.5)
}
P = Delta %*% Q %*% Delta
Sigma <- Lambda %*% P %*% Lambda

Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE,doSym = TRUE)$mat)

nu = n.counties + 1
Q2 = rinvwishart(nu, I)    # note output is an array
Delta2 = matrix(0, n.counties, n.counties)

for (j in 1:n.counties){
  Delta2[j,j] = Q2[j,j]^(-0.5)
}
P2 = Delta2 %*% Q2 %*% Delta2
Sigma2 = Lambda2 %*% P2 %*% Lambda2

n.hunt.i <- ifelse(is.na(chuk_hunt), floor(mean(as.matrix(chuk_hunt), na.rm = T)), NA)
n.harv.i <- as.data.frame(ifelse(is.na(chuk_harv), floor(mean(as.matrix(chuk_harv), na.rm = T)), NA))

Ni <- matrix(NA, nrow = nrow(n.harv.i), ncol = ncol(n.harv.i))
Ni[,1] <- n.harv.i[,1] + 50

## Chukar Site Abundance
chukar_na <- survey.abun

for(i in 1:nrow(chukar_na)){
  for(j in 1:ncol(chukar_na)){
    if(is.na(chukar_na[i,j])){
      chukar_na[i,j] <- floor(mean(as.matrix(survey.abun[i,]), na.rm = T))
    }else{
      chukar_na[i,j] <- NA
    }
  }
}

r.chuk.init <- matrix(NA, nrow = nrow(survey.abun), ncol = ncol(survey.abun)-1)
for(i in 1:nrow(chukar_na)){
  for(j in 2:ncol(chukar_na)){
    r.chuk.init[i,j-1] <- chukar_na[i,j]/chukar_na[i,j-1]
  }
}

C.chuk.init <- survey.abun
C.chuk.init[is.na(C.chuk.init)] <- floor(mean(as.matrix(survey.abun), na.rm = T))
C.chuk.init[,2:ncol(C.chuk.init)] <- NA

econ.inits <- econ_data %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))
bbs.inits <- as.data.frame(bbs.df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA))
awssi.inits <- as.matrix(as.data.frame(awssi.df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))
pdsi.inits <- as.matrix(as.data.frame(pdsi_df) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))
bobcat.inits <- as.matrix(as.data.frame(bobcat.df[,2]) %>% mutate_all(function(x) ifelse(is.na(x), 0, NA)))[,1]

