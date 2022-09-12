lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr"), require, character.only = T)


load(file = "./Output/NDOW_ChukOnly_SEM_output.rdata")
mcmcList1 <- files[[1]]
county_order <- files[[5]]


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
  geom_ribbon(aes(ymin = Est.H.LCL*1000, ymax = Est.H.UCL*1000), alpha = .4, linetype = "dashed") +
  geom_line(aes(y = Est.H*1000, color = County), size = 1) +
  geom_point(aes(y = Obs.H, color = County), size = 1.5) +
  facet_wrap(vars(County), scales = "free_y", ncol = 5) +
  theme_classic(base_size = 18) + 
  labs(y = "Hunter Participation (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10') +
  # scale_color_manual(values = colors.graph) +
  # scale_fill_manual(values = colors.graph) +
  theme(legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(graph.H,  filename = './Graphs/CHUK - Hunter Participation Graph.jpeg',
       dpi = 300, width = 20, height = 15)


graph.N <- ggplot(data = Est.df.c, aes(x = Year, fill = County, color = County)) +
  geom_ribbon(aes(ymin = Est.N.LCL*1000, ymax = Est.N.UCL*1000), alpha = .4, linetype = "dashed") +
  geom_line(aes(y = Est.N*1000, color = County), size = 1) +
  geom_point(aes(y = Obs.N, color = County), size = 1.5) +
  facet_wrap(vars(County), scales = "free_y", ncol = 5) +
  theme_classic(base_size = 18) + 
  labs(y = "Total Harvest (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10', labels = comma) +
  # scale_color_manual(values = colors.graph) +
  # scale_fill_manual(values = colors.graph) +
  theme(legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(graph.N,  filename = './Graphs/CHUK - Total Harvest Graph.jpeg',
       dpi = 300, width = 20, height = 15)


graph.BPH <- ggplot(data = Est.df.c, aes(x = Year, fill = County, color = County)) +
  geom_ribbon(aes(ymin = Est.BPH.LCL, ymax = Est.BPH.UCL), alpha = .4, linetype = "dashed") +
  geom_line(aes(y = Est.BPH, color = County), size = 1) +
  geom_point(aes(y = Obs.BPH, color = County), size = 1.5) +
  facet_wrap(vars(County), scales = "free_y", ncol = 5) +
  theme_classic(base_size = 18) + 
  labs(y = "Birds Harvested per Hunter (log scaled)") +
  scale_x_continuous(breaks= pretty_breaks(n = 4))  + 
  scale_y_continuous(trans='log10', labels = comma) +
  # scale_color_manual(values = colors.graph) +
  # scale_fill_manual(values = colors.graph) +
  theme(legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(graph.BPH,  filename = './Graphs/CHUK - Birds per Hunter Graph.jpeg',
       dpi = 300, width = 20, height = 15)

# harvest correlation
rho.harv.est <- MCMCsummary(mcmcList1, 'rho.harv') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.harv'))) %>%
  mutate(County1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         County2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County1, County2, Estimate = mean) %>%
  pivot_wider(names_from = "County2", values_from = "Estimate") %>%
  dplyr::select(-County1)

rho.harv.plot <- ggcorrplot::ggcorrplot(rho.harv.est, lab = T) +
  labs(title = "Harvest Correlation")
ggsave(rho.harv.plot, filename = './Output/CheckPlot - Rho Harv ChukOnly.jpg',
       dpi = 300, width = 10, height = 10)

rho.hunt.est <- MCMCsummary(mcmcList1, 'rho.hunt') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'rho.hunt'))) %>%
  mutate(County1 = as.numeric(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         County2 = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County1, County2, Estimate = mean) %>%
  pivot_wider(names_from = "County2", values_from = "Estimate") %>%
  dplyr::select(-County1)

rho.hunt.plot <- ggcorrplot::ggcorrplot(rho.hunt.est, lab = T) +
  labs(title = "Hunter Effort Correlation")
ggsave(rho.hunt.plot, filename = './Output/CheckPlot - Rho Hunt ChukOnly.jpg',
       dpi = 300, width = 10, height = 10)