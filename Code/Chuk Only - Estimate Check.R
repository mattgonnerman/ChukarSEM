lapply(c("dplyr", "ggplot2", "coda", "MCMCvis", "stringr"), require, character.only = T)


# load(file = "./Output/NDOW_ChukOnly_SEM_output.rdata")
# mcmcList1 <- files[[1]]
# county_order <- files[[5]]

# Extract important values and plot
test.bph <- MCMCsummary(mcmcList1, 'BPH') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'BPH'))) %>%
  mutate(County = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')
test.H   <- MCMCsummary(mcmcList1, 'H') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'H'))) %>%
  mutate(County = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')
test.N   <- MCMCsummary(mcmcList1, 'N') %>%
  mutate(RowID = rownames(MCMCsummary(mcmcList1, 'N'))) %>%
  mutate(County = as.factor(str_extract(RowID, "(?<=\\[).*?(?=\\,)")),
         Year = as.numeric(str_extract(RowID, "(?<=\\, ).*?(?=\\])"))) %>%
  dplyr::select(County, Year, Estimate = '50%', LCL = '2.5%', UCL = '97.5%')
write.csv(test.bph, file = "./Output/CheckPlot - out_BPH_chukonly.csv", row.names = F)
write.csv(test.H, file = "./Output/CheckPlot - out_H_chukonly.csv", row.names = F)
write.csv(test.N, file =  "./Output/CheckPlot - out_N_chukonly.csv", row.names = F)

est.N <- test.N %>%
  mutate(Year = Year + 1989,
         County = county_order[as.numeric(as.character(County))]) %>%
  dplyr::rename(Est.N = Estimate, Est.N.LCL = LCL, Est.N.UCL = UCL) %>%
  filter(!is.na(County))

est.H <- test.H %>%
  mutate(Year = Year + 1989,
         County = county_order[as.numeric(as.character(County))]) %>%
  dplyr::rename(Est.H = Estimate, Est.H.LCL = LCL, Est.H.UCL = UCL) %>%
  filter(!is.na(County))

est.BPH <- test.bph %>%
  mutate(Year = Year + 1989,
         County = county_order[as.numeric(as.character(County))]) %>%
  dplyr::rename(Est.BPH = Estimate, Est.BPH.LCL = LCL, Est.BPH.UCL = UCL) %>%
  filter(!is.na(County))  

est.check.df <- merge(est.N, est.H, by = c("County", "Year"), all = T)
est.check.df <- merge(est.check.df, est.BPH, by = c("County", "Year"), all = T)

chukharv_data <- read.csv('./Data/ChukarHarvestData.csv') %>%
  mutate(N = round(N),
         H = round(H)) %>%
  filter(Year > cutoff.y) %>%
  arrange(County) %>%
  filter(County %in% unique(est.check.df$County)) %>%
  dplyr::rename(N.obs = N, H.obs = H) %>%
  mutate(BPH.obs = N.obs/H.obs)

est.check.df <- merge(est.check.df, chukharv_data, by = c("County", "Year"), all = T) %>%
  arrange(County, Year)
pred.only.df <- est.check.df %>% filter(Year < 2018 & Year >= cutoff.y) %>%
  mutate(Year = as.factor(Year))

est.check.N <- ggplot(data = pred.only.df, aes(x = Year, group = County)) +
  geom_errorbar(aes(ymin = Est.N.LCL*1, ymax = Est.N.UCL*1, color = County),
                width = .3, size = 1,
                position = position_dodge(width = .5), show.legend = F) +
  geom_point(aes(y = Est.N*1, color = County), shape = 19, size = 4,
             position = position_dodge(width = .5)) +
  geom_point(aes(y = N.obs, group = County), color = "black", shape = 17, size = 3,
             position = position_dodge(width = .5), show.legend = F) +
  theme_classic(base_size = 18) +
  labs(title = "Total Harvest" ) +
  theme(axis.title = element_blank(), legend.position = "right")

est.check.H <- ggplot(data = pred.only.df, aes(x = Year, group = County)) +
  geom_errorbar(aes(ymin = Est.H.LCL*1, ymax = Est.H.UCL*1, color = County),
                width = .3, size = 1,
                position = position_dodge(width = .5), show.legend = F) +
  geom_point(aes(y = Est.H*1, color = County), shape = 19, size = 4,
             position = position_dodge(width = .5)) +
  geom_point(aes(y = H.obs, group = County), color = "black", shape = 17, size = 3,
             position = position_dodge(width = .5), show.legend = F) +
  theme_classic(base_size = 18) +
  labs(title = "Hunter Effort" ) +
  theme(axis.title = element_blank(), legend.position = "right")

est.check.BPH <- ggplot(data = pred.only.df, aes(x = Year, group = County)) +
  geom_errorbar(aes(ymin = Est.BPH.LCL, ymax = Est.BPH.UCL, color = County),
                width = .3, size = 1,
                position = position_dodge(width = .5), show.legend = F) +
  geom_point(aes(y = Est.BPH, color = County), shape = 19, size = 4,
             position = position_dodge(width = .5)) +
  geom_point(aes(y = BPH.obs, group = County), color = "black", shape = 17, size = 3,
             position = position_dodge(width = .5), show.legend = F) +
  theme_classic(base_size = 18) +
  labs(title = "Birds per Hunter" ) +
  theme(axis.title = element_blank(), legend.position = "right")

ggsave(est.check.N, filename = './Output/EstCheck - ChukOnly - N.jpeg',
       dpi = 300, width = 10 + (year.hold - cutoff.y), height = 8)

ggsave(est.check.H, filename = './Output/EstCheck - ChukOnly - H.jpeg',
       dpi = 300, width = 10 + (year.hold - cutoff.y), height = 8)

ggsave(est.check.BPH, filename = './Output/EstCheck - ChukOnly - BPH.jpeg',
       dpi = 300, width = 10 + (year.hold - cutoff.y), height = 8)

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

