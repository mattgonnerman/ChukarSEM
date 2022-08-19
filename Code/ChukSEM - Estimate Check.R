lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr"), require, character.only = T)

load(file = "./Output/NDOW_Upland_SEM_output.rdata")
mcmcList1 <- files[[1]]
mcmcList2 <- files[[1]]

# Extract important values and plot
test.bph <- MCMCsummary(mcmcList1, 'BPH') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'BPH'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')
test.H   <- MCMCsummary(mcmcList1, 'H') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'H'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')
test.N   <- MCMCsummary(mcmcList1, 'N') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'N'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')
write.csv(test.bph, file = "./Output/CheckPlot - out_BPH_all.csv", row.names = F)
write.csv(test.H, file = "./Output/CheckPlot - out_H_all.csv", row.names = F)
write.csv(test.N, file =  "./Output/CheckPlot - out_N_all.csv", row.names = F)

checkplotN <- ggplot(data = test.N, aes(x = Year, y = Estimate, group = as.factor(Species))) +
  geom_line(aes(color = as.factor(Species))) +
  facet_wrap(vars(Region)) +
  scale_y_continuous(trans = "log10") +
  labs(title = "Total Harvest")
ggsave(checkplotN, filename = "./Output/CheckPlot - N.jpg", dpi = 300)

checkplotH <- ggplot(data = test.H, aes(x = Year, y = Estimate, as.factor(Species))) +
  geom_line(aes(color = as.factor(Species))) +
  facet_wrap(vars(Region)) +
  scale_y_continuous(trans = "log10") +
  labs(title = "Hunter Effort")
ggsave(checkplotH, filename = "./Output/CheckPlot - H.jpg", dpi = 300)

checkplotBPH <- ggplot(data = test.bph, aes(x = Year, y = Estimate, as.factor(Species))) +
  geom_line(aes(color = as.factor(Species))) +
  facet_wrap(vars(Region)) +
  scale_y_continuous(trans = "log10") +
  labs(title = "Birds Per Hunter")
ggsave(checkplotBPH, filename = "./Output/CheckPlot - BPH.jpg", dpi = 300)


### Check Model Outputs against observed values

check.species <- species

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
  filter(Year <= final.y) %>%
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
pred.only.df <- est.check.df %>% filter(Year >= cutoff.y) %>%
  mutate(Year = as.factor(Year))

write.csv(as.data.frame(pred.only.df), "./Output/EstCheck - Estimates_N_H_BPH.csv", row.names = F)

est.check.N <- ggplot(data = pred.only.df, aes(x = Year, group = Region)) +
  geom_errorbar(aes(ymin = Est.N.LCL*100, ymax = Est.N.UCL*100, color = Region),
                width = .3, size = 1,
                position = position_dodge(width = .5), show.legend = F) +
  geom_point(aes(y = Est.N*100, color = Region), shape = 19, size = 4,
             position = position_dodge(width = .5)) +
  geom_point(aes(y = Obs.N, group = Region), color = "black", shape = 17, size = 3,
             position = position_dodge(width = .5), show.legend = F) +
  facet_wrap(vars(Species), scales = "free_y", nrow = 2) +
  theme_classic(base_size = 18) +
  labs(title = "Total Harvest" ) +
  scale_color_manual(name = "Region", values = c("#0033a0", "#ffa300"), labels = c("Eastern", "Western")) +
  theme(axis.title = element_blank(), legend.position = c(.9,.25))

ggsave(est.check.N, filename = './Output/EstCheck - N.jpeg',
       dpi = 300, width = 10 + (year.hold - cutoff.y), height = 8)


est.check.H <- ggplot(data = pred.only.df, aes(x = Year, group = Region)) +
  geom_errorbar(aes(ymin = Est.H.LCL*100, ymax = Est.H.UCL*100, color = Region),
                width = .3, size = 1,
                position = position_dodge(width = .5), show.legend = F) +
  geom_point(aes(y = Est.H*100, color = Region), shape = 19, size = 4,
             position = position_dodge(width = .5)) +
  geom_point(aes(y = Obs.H, group = Region), color = "black", shape = 17, size = 3,
             position = position_dodge(width = .5), show.legend = F) +
  facet_wrap(vars(Species), scales = "free_y", nrow = 2) +
  theme_classic(base_size = 18) +
  labs(title = "Hunter Participation" ) +
  scale_color_manual(name = "Region", values = c("#0033a0", "#ffa300"), labels = c("Eastern", "Western")) +
  theme(axis.title = element_blank(), legend.position = c(.9,.25))

ggsave(est.check.H, filename = './Output/EstCheck - H.jpeg',
       dpi = 300, width = 10 + (year.hold - cutoff.y), height = 8)


est.check.BPH <- ggplot(data = pred.only.df, aes(x = Year, group = Region)) +
  geom_errorbar(aes(ymin = Est.BPH.LCL, ymax = Est.BPH.UCL, color = Region),
                width = .3, size = 1,
                position = position_dodge(width = .5), show.legend = F) +
  geom_point(aes(y = Est.BPH, color = Region), shape = 19, size = 4,
             position = position_dodge(width = .5)) +
  geom_point(aes(y = Obs.BPH, group = Region), color = "black", shape = 17, size = 3,
             position = position_dodge(width = .5), show.legend = F) +
  facet_wrap(vars(Species), scales = "free_y", nrow = 2) +
  theme_classic(base_size = 18) +
  labs(title = "Birds per Hunters" ) +
  scale_color_manual(name = "Region", values = c("#0033a0", "#ffa300"), labels = c("Eastern", "Western")) +
  theme(axis.title = element_blank(), legend.position = c(.9,.25))

ggsave(est.check.BPH, filename = './Output/EstCheck - BPH.jpeg',
       dpi = 300, width = 10 + (year.hold - cutoff.y), height = 8)



# harvest correlation
rho.harv.est <- MCMCsummary(mcmcList1, 'rho.harv') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.harv'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 7, ncol = 7))
rho.harv.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.harv.est)){
  rho.harv.list[[rho.harv.est$Region[i]]][rho.harv.est$Species1[i], rho.harv.est$Species2[i]] <- rho.harv.est$Estimate[i]
}

rho.harv.e <- rho.harv.list[[1]]
rho.harv.w <- rho.harv.list[[2]]
write.csv(rho.harv.e, file = "./Output/CheckPlot - rho_harv_east.csv")
write.csv(rho.harv.w, file = "./Output/CheckPlot - rho_harv_west.csv")

rho.harv.e.plot <- ggcorrplot::ggcorrplot(rho.harv.e, lab = T) +
  labs(title = "Harvest Correlation - East")
ggsave(rho.harv.e.plot, filename = './Output/CheckPlot - Rho Harv East.jpg', dpi = 300)
rho.harv.w.plot <- ggcorrplot::ggcorrplot(rho.harv.w, lab = T) +
  labs(title = "Harvest Correlation - West")
ggsave(rho.harv.w.plot, filename = './Output/CheckPlot - Rho Harv West.jpg', dpi = 300)

# hunter correlation
rho.hunt.est <- MCMCsummary(mcmcList1, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 7, ncol = 7))
rho.hunt.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.hunt.est)){
  rho.hunt.list[[rho.hunt.est$Region[i]]][rho.hunt.est$Species1[i], rho.hunt.est$Species2[i]] <- rho.hunt.est$Estimate[i]
}

rho.hunt.e <- rho.hunt.list[[1]]
rho.hunt.w <- rho.hunt.list[[2]]
write.csv(rho.hunt.e, file = "./Output/CheckPlot - rho_hunt_east.csv")
write.csv(rho.hunt.w, file = "./Output/CheckPlot - rho_hunt_west.csv")

rho.hunt.e.plot <- ggcorrplot::ggcorrplot(rho.hunt.e, lab = T) +
  labs(title = "Hunter Correlation - East")
ggsave(rho.hunt.e.plot, filename = './Output/CheckPlot - Rho Hunt East.jpg', dpi = 300)
rho.hunt.w.plot <- ggcorrplot::ggcorrplot(rho.hunt.w, lab = T) +
  labs(title = "Hunter Correlation - West")
ggsave(rho.hunt.w.plot, filename = './Output/CheckPlot - Rho Hunt West.jpg', dpi = 300)




### Compare Beta Coefficients among species
# beta.co.vec <- ()
# 
# for(i in 1:length(beta.co.vec)){
#   rho.harv.est <- MCMCsummary(mcmcList1, beta.co.vec[i]) %>%
#     mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.harv'))) %>%
#     mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
#            Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
#            Region = sub('.*\\,', '', RowID)) %>%
#     mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
#     dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')
#  
# }
 