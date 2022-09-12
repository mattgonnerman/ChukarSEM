lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr"), require, character.only = T)


load(file = "./Output/NDOW_ChukOnly_SEM_output.rdata")
mcmcList1 <- files[[1]]
mcmcList2 <- files[[2]]
county_order <- files[[5]]
color_county <- c(viridis::viridis(n =13), "#404040")
county_list <- c(county_order, "Statewide")
colors.graph <- setNames(color_county,  county_coef)

# Extract important values and plot
Est.BPH.c <- MCMCsummary(mcmcList1, 'BPH') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'BPH'))) %>%
  mutate(County = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Estimate = '50%', LCL = '2.5%', UCL = '97.5%') %>%
  mutate(Year = Year + 1989,
         County = county_order[as.numeric(as.character(County))]) %>%
  dplyr::rename(Est.BPH = Estimate, Est.BPH.LCL = LCL, Est.BPH.UCL = UCL) %>%
  filter(!is.na(County))  


Est.H.c <- MCMCsummary(mcmcList1, 'H') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'H'))) %>%
  mutate(County = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')%>%
  mutate(Year = Year + 1989,
         County = county_order[as.numeric(as.character(County))]) %>%
  dplyr::rename(Est.H = Estimate, Est.H.LCL = LCL, Est.H.UCL = UCL) %>%
  filter(!is.na(County))


Est.N.c <- MCMCsummary(mcmcList1, 'N') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'N'))) %>%
  mutate(County = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Estimate = '50%', LCL = '2.5%', UCL = '97.5%') %>%
  mutate(Year = Year + 1989,
         County = county_order[as.numeric(as.character(County))]) %>%
  dplyr::rename(Est.N = Estimate, Est.N.LCL = LCL, Est.N.UCL = UCL) %>%
  filter(!is.na(County))

chukharv_data <- read.csv('./Data/ChukarHarvestData.csv') %>%
  mutate(N = round(N),
         H = round(H)) %>%
  filter(Year %in% Est.N.c$Year) %>%
  arrange(County) %>%
  filter(County %in% county_order) %>%
  dplyr::rename(Obs.N = N, Obs.H = H) %>%
  mutate(Obs.BPH = Obs.N/Obs.H)

Est.df.c <- merge(Est.N.c, Est.H.c, by = c("County", "Year"), all = T)
Est.df.c <- merge(Est.df.c, Est.BPH.c, by = c("County", "Year"), all = T)
Est.df.c <- merge(Est.df.c, chukharv_data, by = c("County", "Year"), all = T) %>%
  arrange(County, Year)


graph.H <- ggplot(data = Est.df.c, aes(x = Year, fill = County, color = County)) +
  geom_ribbon(aes(ymin = Est.H.LCL*1000, ymax = Est.H.UCL*1000), alpha = .3, linetype = "dashed") +
  geom_line(aes(y = Est.H*1000, color = County), size = 1) +
  geom_point(aes(y = Obs.H), color = "black", size = 1.5, show.legend = F) +
  facet_wrap(vars(County), scales = "free_y", ncol = 5) +
  theme_classic(base_size = 18) + 
  labs(y = "Hunter Participation (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = colors.graph[1:13]) +
  scale_fill_manual(values = colors.graph[1:13]) +
  theme(legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(graph.H,  filename = './Graphs/CHUK - Hunter Participation Graph.jpeg',
       dpi = 300, width = 20, height = 15)


graph.N <- ggplot(data = Est.df.c, aes(x = Year, fill = County, color = County)) +
  geom_ribbon(aes(ymin = Est.N.LCL*1000, ymax = Est.N.UCL*1000), alpha = .3, linetype = "dashed") +
  geom_line(aes(y = Est.N*1000, color = County), size = 1) +
  geom_point(aes(y = Obs.N), color = "black", size = 1.5, show.legend = F) +
  facet_wrap(vars(County), scales = "free_y", ncol = 5) +
  theme_classic(base_size = 18) + 
  labs(y = "Total Harvest (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10', labels = comma) +
  scale_color_manual(values = colors.graph[1:13]) +
  scale_fill_manual(values = colors.graph[1:13]) +
  theme(legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(graph.N,  filename = './Graphs/CHUK - Total Harvest Graph.jpeg',
       dpi = 300, width = 20, height = 15)


graph.BPH <- ggplot(data = Est.df.c, aes(x = Year, fill = County, color = County)) +
  geom_ribbon(aes(ymin = Est.BPH.LCL, ymax = Est.BPH.UCL), alpha = .3, linetype = "dashed") +
  geom_line(aes(y = Est.BPH, color = County), size = 1) +
  geom_point(aes(y = Obs.BPH), color = "black", size = 1.5, show.legend = F) +
  facet_wrap(vars(County), scales = "free_y", ncol = 5) +
  theme_classic(base_size = 18) + 
  labs(y = "Birds Harvested per Hunter (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10', labels = comma) +
  scale_color_manual(values = colors.graph[1:13]) +
  scale_fill_manual(values = colors.graph[1:13]) +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(graph.BPH,  filename = './Graphs/CHUK - Birds per Hunter Graph.jpeg',
       dpi = 300, width = 20, height = 15)



### Correlation
rho.harv.est <- MCMCsummary(mcmcList1, 'rho.harv') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.harv'))) %>%
  mutate(County1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         County2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County1, County2, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

rho.harv <- as.data.frame(matrix(NA, nrow = 13, ncol = 13))
for(i in 1:nrow(rho.harv.est)){
  rho.harv[rho.harv.est$County1[i], rho.harv.est$County2[i]] <- rho.harv.est$Estimate[i]
}
colnames(rho.harv) <- county_order
rownames(rho.harv) <- county_order

rho.hunt.est <- MCMCsummary(mcmcList1, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.hunt'))) %>%
  mutate(County1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         County2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County1, County2, Estimate = mean, LCL = '2.5%', UCL = '97.5%')

rho.hunt <- as.data.frame(matrix(NA, nrow = 13, ncol = 13))
for(i in 1:nrow(rho.hunt.est)){
  rho.hunt[rho.hunt.est$County1[i], rho.hunt.est$County2[i]] <- rho.hunt.est$Estimate[i]
}
colnames(rho.hunt) <- county_order
rownames(rho.hunt) <- county_order

rho.harv.plot <- ggcorrplot::ggcorrplot(rho.harv, lab = T, type = "upper",
                                          ggtheme = ggplot2::theme_classic,
                                          colors = c("#6D9EC1", "white", "#E46726")) +
  labs(title = "Total Harvest") +
  theme(legend.position = "none")
rho.hunt.plot <- ggcorrplot::ggcorrplot(rho.hunt, lab = T, type = "upper",
                                          ggtheme = ggplot2::theme_classic,
                                          colors = c("#6D9EC1", "white", "#E46726")) +
  labs(title = "Hunter Participation") +
  theme(legend.position = "none")

require(patchwork)
rho.plot <- rho.hunt.plot / rho.harv.plot

ggsave(rho.plot, filename = './Graphs/CHUK - Correlation Plots.jpg',
       dpi = 300, width = 10, height = 15)


### Beta Coefficients
reg.coef <- MCMCsummary(mcmcList2, 'beta.econ.hunt') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "Economy") %>%
  dplyr::select(County, ID, mean, sd) %>%
  mutate(County = as.factor(county_order[County]))

reg.coef <- MCMCsummary(mcmcList2, 'beta.wintsev.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = "Statewide",
         ID = "Winter") %>%
  dplyr::select(County, ID, mean, sd) %>%
  rbind(reg.coef, .)

reg.coef <- MCMCsummary(mcmcList2, 'beta.bbs.harv')%>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = "Statewide",
         ID = "Predator") %>%
  dplyr::select(County, ID, mean, sd) %>%
  rbind(reg.coef, .)

reg.coef <- MCMCsummary(mcmcList2, 'beta.hunter.harv') %>%
  mutate(RowID = rownames(.)) %>% `rownames<-`( NULL ) %>%
  mutate(County = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\])")),
         ID = "Hunter") %>%
  dplyr::select(County, ID, mean, sd) %>%
  mutate(County = as.factor(county_order[County])) %>%
  rbind(reg.coef, .) %>%
  rename(Estimate = mean)
reg.coef$County <- reorder(reg.coef$County)

levels(reg.coef$County)

chuk.coef.plot <- ggplot(reg.coef, aes(y = Estimate, x = as.factor(ID))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.5) +
  geom_errorbar(aes(ymin = Estimate - 1.96*sd,
                    ymax = Estimate + 1.96*sd,
                    color = County), size = 1, width = .25, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = Estimate,
                 color = County), size = 2, position = position_dodge(width = 0.5)) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = colors.graph) +
  labs(y = "Regression Coefficient Value", x = element_blank()) +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))

ggsave(chuk.coef.plot, filename = './Graphs/CHUK - Coefficient Plots.jpg',
       dpi = 300, width = 8, height = 8)
