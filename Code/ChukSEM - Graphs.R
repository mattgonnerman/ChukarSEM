lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr", "scales"), require, character.only = T)


### Load Estimated Values/CI
load(file = "./Output/NDOW_Upland_SEM_output.rdata")
mcmcList1 <- files[[1]]
mcmcList2 <- files[[2]]
species.constant <- species <- check.species <- sort(c("SAGR", "CHUK", "BLGR", "CAQU", "HUPA", "PHEA"))
colors.constant <- c("#E8A323","#23E8A3","#A323E8","#074A36","#36074A", "#4A3607")
colors.graph <- setNames(colors.constant,  species.constant)


est.BPH <- MCMCsummary(mcmcList1, 'BPH') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'BPH'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%') %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.BPH = Estimate, Est.BPH.LCL = LCL, Est.BPH.UCL = UCL) %>%
  filter(!is.na(Species))  
est.BPH2 <- MCMCsummary(mcmcList1, 'BPH2') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'BPH2'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%') %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.BPH = Estimate, Est.BPH.LCL = LCL, Est.BPH.UCL = UCL) %>%
  filter(!is.na(Species))  
est.H   <- MCMCsummary(mcmcList1, 'H') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'H'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%') %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.H = Estimate, Est.H.LCL = LCL, Est.H.UCL = UCL) %>%
  filter(!is.na(Species))
est.N   <- MCMCsummary(mcmcList1, 'N') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'N'))) %>%
  mutate(Species = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species, Year, Region, Estimate = '50%', LCL = '2.5%', UCL = '97.5%') %>%
  mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
         Year = Year + 1975,
         Species = check.species[as.numeric(Species)]) %>%
  dplyr::rename(Est.N = Estimate, Est.N.LCL = LCL, Est.N.UCL = UCL) %>%
  filter(!is.na(Species))

### Load Observed Values
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

est.df <- merge(est.N, est.H, by = c("Region", "Species", "Year"), all = T)
est.df <- merge(est.df, est.BPH, by = c("Region", "Species", "Year"), all = T)
est.df <- merge(est.df, obs.HarvestData, by = c("Region", "Species", "Year"), all = T) %>%
  arrange(Species, Year, Region)

graph.H <- ggplot(data = est.df, aes(x = Year, fill = Species, color = Species)) +
  geom_ribbon(aes(ymin = Est.H.LCL*1000, ymax = Est.H.UCL*1000), alpha = .4, linetype = "dashed") +
  geom_line(aes(y = Est.H*1000, color = Species), size = 1) +
  geom_point(aes(y = Obs.H, color = Species), size = 1.5) +
  facet_wrap(vars(Region), scales = "free_y", nrow = 2) +
  theme_classic(base_size = 18) + 
  labs(y = "Hunter Participation (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = colors.graph) +
  scale_fill_manual(values = colors.graph) +
  theme(legend.position = "right")+
  guides(color=guide_legend(nrow=6, byrow=TRUE),
         fill=guide_legend(nrow=6, byrow=TRUE))
ggsave(graph.H,  filename = './Graphs/Full - Hunter Participation Graph.jpeg',
       dpi = 300, width = 11, height = 8.5)


graph.N <- ggplot(data = est.df, aes(x = Year, fill = Species, color = Species)) +
  geom_ribbon(aes(ymin = Est.N.LCL*1000, ymax = Est.N.UCL*1000), alpha = .4, linetype = "dashed") +
  geom_line(aes(y = Est.N*1000, color = Species), size = 1) +
  geom_point(aes(y = Obs.N, color = Species), size = 1.5) +
  facet_wrap(vars(Region), scales = "free_y", nrow = 2) +
  theme_classic(base_size = 18) + 
  labs(y = "Total Harvest (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10',labels = comma) +
  scale_color_manual(values = colors.graph) +
  scale_fill_manual(values = colors.graph) +
  theme(legend.position = "right")+
  guides(color=guide_legend(nrow=6, byrow=TRUE),
         fill=guide_legend(nrow=6, byrow=TRUE))
ggsave(graph.N,  filename = './Graphs/Full - Total Harvest Graph.jpeg',
       dpi = 300, width = 11, height = 8.5)


graph.BPH <- ggplot(data = est.df, aes(x = Year, fill = Region, color = Region)) +
  geom_ribbon(aes(ymin = Est.BPH.LCL, ymax = Est.BPH.UCL, fill = Species), alpha = .4, linetype = "dashed") +
  geom_line(aes(y = Est.BPH, color = Species), size = 1) +
  geom_point(aes(y = Obs.BPH), color = "black", size = 1.5) +
  facet_wrap(~ Region + Species, scales = "free_y", nrow = 4) +
  theme_classic(base_size = 18) + 
  labs(y = "Birds Harvested per Hunter (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4)) + 
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = colors.graph) +
  scale_fill_manual(values = colors.graph) +
  theme(legend.position = "right") +
  guides(color=guide_legend(nrow=6, byrow=TRUE),
         fill=guide_legend(nrow=6, byrow=TRUE))
ggsave(graph.BPH,  filename = './Graphs/Full - Birds per Hunter Graph.jpeg',
       dpi = 400, width = 20, height = 15)


### Correlation
rho.harv.est <- MCMCsummary(mcmcList1, 'rho.harv') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.harv'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

blank.cor.df <- as.data.frame(matrix(NA, nrow = 6, ncol = 6))
rho.harv.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.harv.est)){
  rho.harv.list[[rho.harv.est$Region[i]]][rho.harv.est$Species1[i], rho.harv.est$Species2[i]] <- rho.harv.est$Estimate[i]
}

rho.harv.est %>% filter( Species1 != Species2) %>% filter(Species1 > Species2) %>%
  summarize(Mean = mean(Estimate))

rho.harv.e <- rho.harv.list[[1]]
colnames(rho.harv.e) <- species
rownames(rho.harv.e) <- species
rho.harv.w <- rho.harv.list[[2]]
colnames(rho.harv.w) <- species
rownames(rho.harv.w) <- species

rho.harv.e.plot <- ggcorrplot::ggcorrplot(rho.harv.e, lab = T, type = "upper",
                                          ggtheme = ggplot2::theme_classic,
                                          colors = c("#6D9EC1", "white", "#E46726")) +
  labs(title = "Total Harvest - Eastern") +
  theme(legend.position = "none")
rho.harv.w.plot <- ggcorrplot::ggcorrplot(rho.harv.w, lab = T, type = "upper",
                                          ggtheme = ggplot2::theme_classic,
                                          colors = c("#6D9EC1", "white", "#E46726")) +
  labs(title = "Total Harvest - Western") +
  theme(legend.position = "none")

rho.harv.plot <- rho.harv.e.plot + plot_spacer() + rho.harv.w.plot + plot_layout(guides = "collect", widths = c(1, .1 ,1))



rho.hunt.est <- MCMCsummary(mcmcList1, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.hunt'))) %>%
  mutate(Species1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Species2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.factor(str_sub(Region,1,nchar(Region)-1))) %>%
  dplyr::select(Species1, Species2, Region, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

rho.hunt.est %>% filter( Species1 != Species2) %>% filter(Species1 > Species2) %>%
  summarize(Mean = mean(Estimate))

blank.cor.df <- as.data.frame(matrix(NA, nrow = 6, ncol = 6))
rho.hunt.list <- list(blank.cor.df, blank.cor.df)
for(i in 1:nrow(rho.hunt.est)){
  rho.hunt.list[[rho.hunt.est$Region[i]]][rho.hunt.est$Species1[i], rho.hunt.est$Species2[i]] <- rho.hunt.est$Estimate[i]
}

rho.hunt.e <- rho.hunt.list[[1]]
colnames(rho.hunt.e) <- species
rownames(rho.hunt.e) <- species
rho.hunt.w <- rho.hunt.list[[2]]
colnames(rho.hunt.w) <- species
rownames(rho.hunt.w) <- species

rho.hunt.e.plot <- ggcorrplot::ggcorrplot(rho.hunt.e, lab = T, type = "upper",
                                          ggtheme = ggplot2::theme_classic,
                                          colors = c("#6D9EC1", "white", "#E46726")) +
  labs(title = "Hunter Participation - Eastern") +
  theme(legend.position = "none")
rho.hunt.w.plot <- ggcorrplot::ggcorrplot(rho.hunt.w, lab = T, type = "upper",
                                          ggtheme = ggplot2::theme_classic,
                                          colors = c("#6D9EC1", "white", "#E46726")) +
  labs(title = "Hunter Participation - Western") +
  theme(legend.position = "none")

require(patchwork)
rho.hunt.plot <- rho.hunt.e.plot + plot_spacer() + rho.hunt.w.plot + plot_layout(guides = "collect", widths = c(1, .1 ,1))

rho.plot <- rho.hunt.plot / rho.harv.plot

ggsave(rho.plot, filename = './Graphs/FULL - Correlation Plots.jpg',
       dpi = 300, width = 10, height = 10)


### Beta Coefficients
Hunt.reg.coef1 <- MCMCsummary(mcmcList2, 'alpha.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID))%>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "b.0") %>%
  dplyr::select(Species, Region, A.hunt = mean, A.hunt.sd = sd) 
Hunt.reg.coef <- MCMCsummary(mcmcList2, 'beta.econ.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.ECON") %>%
  dplyr::select(Species, B.econ.hunt = mean, B.econ.hunt.sd = sd) %>%
  merge(., expand.grid(Species = 1:6, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Hunt.reg.coef1, ., by = c("Species", "Region"), all =T)

### Total Harvest Regression Coefficients
Harv.reg.coef1 <- MCMCsummary(mcmcList2, 'alpha.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = sub('.*\\,', '', RowID)) %>%
  mutate(Region = as.numeric(str_sub(Region,1,nchar(Region)-1)),
         ID = "alpha.harv") %>%
  dplyr::select(Species, Region, A.harv = mean, A.harv.sd = sd) 

Harv.reg.coef2 <- MCMCsummary(mcmcList2, 'beta.wintsev.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.WINT") %>%
  dplyr::select(Species, B.ws.harv = mean, B.ws.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:6, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef1, ., by = c("Species", "Region"), all =T)

Harv.reg.coef3 <- MCMCsummary(mcmcList2, 'beta.bbs.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "b.BBS") %>%
  dplyr::select(Species, B.bbs.harv = mean, B.bbs.harv.sd = sd) %>%
  merge(., expand.grid(Species = 1:6, Region = 1:2), by = "Species", all.x = T) %>%
  merge(Harv.reg.coef2, ., by = c("Species", "Region"), all =T)

Harv.reg.coef <- MCMCsummary(mcmcList2, 'beta.hunter.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(Species = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Region = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(Species, Region, B.hunter = mean, B.hunter.sd = sd) %>%
  merge(Harv.reg.coef3, ., by = c("Species", "Region"), all =T)

reg.coef <- merge(Hunt.reg.coef, Harv.reg.coef, by = c("Species", "Region"), all = T) %>%
  mutate(Species = as.factor(check.species[Species]))

B.econ <- reg.coef %>%
  dplyr::select(Species, B.econ.hunt, B.econ.hunt.sd) %>% distinct()

econ.plot <- ggplot(B.econ, aes(x = Species)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.5) +
  geom_errorbar(aes(ymin = B.econ.hunt - 1.96*B.econ.hunt.sd,
                    ymax = B.econ.hunt + 1.96*B.econ.hunt.sd,
                    color = Species), size = 1, width = .25) +
  geom_point(aes(y = B.econ.hunt,
                 color = Species), size = 2) +
  theme_classic(base_size = 12) +
  labs(title = "Hunter Response to Economic Conditions",
       y = "Beta Coefficient", x = element_blank()) +
  scale_color_manual(values = colors.graph) +
  theme(legend.position = "none")


B.bbs <- reg.coef %>%
  dplyr::select(Species, B.bbs.harv, B.bbs.harv.sd) %>% distinct()

bbs.plot <- ggplot(B.bbs, aes(x = Species)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.5) +
  geom_errorbar(aes(ymin = B.bbs.harv - 1.96*B.bbs.harv.sd,
                    ymax = B.bbs.harv + 1.96*B.bbs.harv.sd,
                    color = Species), size = 1, width = .25) +
  geom_point(aes(y = B.bbs.harv,
                 color = Species), size = 2) +
  theme_classic(base_size = 12) +
  labs(title = "Change in Total Harvest by Predator Index",
       y = "Beta Coefficient", x = element_blank()) +
  scale_color_manual(values = colors.graph) +
  theme(legend.position = "none")


B.wint <- reg.coef %>%
  dplyr::select(Species, B.ws.harv, B.ws.harv.sd) %>% distinct()

wint.plot <- ggplot(B.wint, aes(x = Species)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.5) +
  geom_errorbar(aes(ymin = B.ws.harv - 1.96*B.ws.harv.sd,
                    ymax = B.ws.harv + 1.96*B.ws.harv.sd,
                    color = Species), size = 1, width = .25) +
  geom_point(aes(y = B.ws.harv,
                 color = Species), size = 2) +
  theme_classic(base_size = 12) +
  labs(title = "Change in Total Harvest by Winter Severity",
       y = "Beta Coefficient", x = element_blank()) +
  scale_color_manual(values = colors.graph) +
  theme(legend.position = "none")

B.hunt <- reg.coef %>%
  dplyr::select(Species, Region, B.hunter, B.hunter.sd) %>% distinct() %>%
  mutate(Region = as.factor(ifelse(Region == 1, "Western", "Eastern")))

hunter.plot <- ggplot(B.hunt, aes(x = Region)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.5) +
  geom_errorbar(aes(ymin = B.hunter - 1.96*B.hunter.sd,
                    ymax = B.hunter + 1.96*B.hunter.sd,
                    color = Species), size = 1, width = .25, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = B.hunter,color = Species),
             position = position_dodge(width = 0.5), size = 2) +
  theme_classic(base_size = 12) +
  labs(title = "Change in Total Harvest by Hunter Participation",
       y = "Beta Coefficient", x = element_blank()) +
  scale_color_manual(values = colors.graph) +
  theme(legend.position = "right")

require(patchwork)

full.coef <- econ.plot + hunter.plot + bbs.plot + wint.plot + plot_layout(nrow = 2)
ggsave(full.coef, filename = './Graphs/FULL - Coefficient Plots.jpg',
       dpi = 300, width = 10, height = 10)
