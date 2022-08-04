load("ChukarSEM_model_output.rdata")

rho1 <- MCMCpstr(mcmcList2,'rho')
rho2 <- MCMCpstr(mcmcList2,'rho2')

rho1a <- rho1$rho[,,1]
rho1b <- rho1$rho[,,2]
rownames(rho1a) <- colnames(rho1a) <- species[-c(4,8)]
rownames(rho1b) <- colnames(rho1b) <- species[-c(4,8)]

rho1a[lower.tri(rho1a)] <- rho1b[lower.tri(rho1b)]

rho2a <- rho2$rho2[,,1]
rho2b <- rho2$rho2[,,2]
rownames(rho2a) <- colnames(rho2a) <- species[-c(4,8)]
rownames(rho2b) <- colnames(rho2b) <- species[-c(4,8)]

rho2a[lower.tri(rho2a)] <- rho2b[lower.tri(rho2b)]

diag(rho1a) <-NA
diag(rho2a) <-NA

library(ggcorrplot)

ggcorrplot(rho1a, outline.col = "white",  show.diag = FALSE, lab= TRUE)
ggcorrplot(rho2a, outline.col = "white",  show.diag = FALSE, lab= TRUE)

summary_one <- MCMCsummary(mcmcList1)
summary_two <- MCMCsummary(mcmcList2)

test.bph <- MCMCsummary(mcmcList2, 'bph')
test.H   <- MCMCsummary(mcmcList2, 'H')
test.N   <- MCMCsummary(mcmcList2, 'N')

test.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                        year = 1976:2021, 
                                        region = c('west','east') ),
                            test.H, hunters = c(abind(hunters, array(NA, dim = c(7,4,2)),along = 2)))

test.dfn <- cbind.data.frame(expand.grid(species = species[-c(4,8)], 
                                         year = 1976:2021, 
                                         region = c('west','east') ),
                             test.N, upland = c(abind(upland, array(NA, dim = c(7,4,2)),along = 2)))

fig_upland <-   ggplot(data = test.dfn, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`),  size = .75) +
  geom_pointrange2(data = subset(test.dfn,year>2017), aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), color = 'red',  size = .75) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1000, color = NA) +  geom_hline(yintercept = 0, color = NA) +
  labs(y = 'Predicted number of birds reported as harvested', x = 'Year') +
  scale_fill_flat_d() + theme_modern()
fig_upland

fig_hunters <-  ggplot(data = test.df, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`),  size = .75) +
  geom_pointrange2(data = subset(test.df,year>2017), aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), color = 'red',  size = .75) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 0) + labs(y = 'Predicted number of hunters', x = 'Year') +
  scale_fill_flat_d() + theme_modern()
fig_hunters

ggsave(fig_hunters, file = 'fig_hunters.tiff', height = 6, width = 15, dpi = 320)

ggplot(data = test.df, aes(x = year, y = mean)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), fill = 'black', shape = 21, size = .75, position = position_dodge(1)) +
  geom_point(aes(x = year, y = hunters), color = 'red') +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1) +
  scale_fill_flat_d() + theme_modern()

test <- MCMCsummary(mcmcList2, 'pred')
test.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], year = 1977:2021, region = c('east','west') ),
                            test, hunters =  c(abind(hunters, array(NA, dim = c(7,3,2)),along = 2)))

fig1 <- ggplot(data = test.df, aes(x = hunters, y = mean, group = region)) +
  geom_pointrange2(aes(x =hunters, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Pressure', x = 'Year') + facet_wrap(region ~ species, scales = 'free', ncol = 7) +  scale_fill_flat_d() + theme_modern(legend.position = 'none')
fig1

test1 <- MCMCsummary(mcmcList2, 'pred1')
test1.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)],  region = c('east','west'), year = 1976:2021 ), test1)

fig2 <- ggplot(data = test1.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Numbers', x = 'Year') + facet_grid(region ~ species) +  scale_fill_flat_d() + theme_modern()
fig2

test.NH <- cbind.data.frame(expand.grid(species = species[-c(4,8)], year = 1976:2021, region = c('west','east')), test.bph)

fig1 <- ggplot(data = subset(test.NH, year == 2017), aes(x = species, y = mean, group = region)) +
  geom_pointrange2(aes(x =species, y = mean, ymin = `2.5%`,ymax = `97.5%`, shape = region),color = 'black', size = .75, position = position_dodge(1)) +
  # geom_point(aes(x = species, y = bph,shape = region),size = 2, color = 'red', position = position_dodge(1)) +
  geom_text(data = subset(test.dfn, year == 2017), aes(x =species, y = -1, group = region, label = round(mean,0)),position = position_dodge(1))+
  geom_text(data = subset(test.dfn, year == 2017), aes(x =species, y = -3, group = region, label = round(upland,0)),position = position_dodge(1), color = 'red')+
  geom_hline(yintercept = 0) + labs(y = 'Forecasted number of individuals harvested per hunter in 2017', x = 'Species') +
  theme_modern()
fig1
ggsave(fig1, file = 'prelim_model.tiff', dpi = 320, height = 6, width = 10)

ggplot(data = test.NH, aes(x = year, y = mean, fill = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, shape = region),color = 'black', size = .75, position = position_dodge(1)) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  geom_hline(yintercept = 1) +
  scale_fill_flat_d() + theme_modern()

ggplot(data = test.NH, aes(x = year, y = mean, fill = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`), shape = 21, size = .75, position = position_dodge(1)) +
  facet_wrap(region ~ species, ncol = 7, scales = 'free_y') +
  geom_hline(yintercept = 0) + labs(y = 'Estimated number of individuals harvested per hunter', x = 'Species') +
  scale_fill_flat_d() +  scale_color_metro_d() + theme_modern()

library(see)
ggplot(data = test.df, aes(x = hunters, y = mean, fill = region)) +
  geom_pointrange2(aes(x = hunters, y = mean, ymin = `2.5%`,ymax = `97.5%`), linetype = 'dashed', shape = 21, size = .66, position = position_dodge(1)) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  labs(x = 'Number of self-reported active hunters in year t', y = 'Observed impact on the number of individuals reportd as harvested in year t+1') +
  geom_hline(yintercept = 0) +
  scale_fill_flat_d() + theme_modern()

ggplot(data = test.df, aes(x = year, y = mean, fill = region)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`,ymax = `97.5%`), linetype = 'dashed', shape = 21, size = .66, position = position_dodge(1)) +
  facet_wrap(region~species,scales = 'free', ncol = 7) +
  labs(x = 'Year', y = 'Observed impact of active hunters in year t on the number\n of individuals reported as harvested in year t+1') +
  geom_hline(yintercept = 0) +
  scale_fill_flat_d() + theme_modern()

plot(hunters[1,1:41,1], test_matrix[1,,1])


N_bird   <- MCMCsummary(mcmcList2, 'N')
N_hunt   <- MCMCsummary(mcmcList2, 'H')
N_survey <- MCMCsummary(mcmcList2, 'C')

N_bird$Species <- rep(species[-c(4,8)], times = 46 * 2)
N_bird$Year    <- rep(1976:2021, each = 7, times = 2)
N_bird$Region  <- rep(c("West",'East'), each = 7 *  46)
N_bird$Type <- rep('Bird')

N_hunt$Species <- rep(species[-c(4,8)], times = 46 * 2)
N_hunt$Year    <- rep(1976:2021, each = 7, times = 2)
N_hunt$Region  <- rep(c("West",'East'), each = 7 *  46)
N_hunt$Type <- rep('Hunter')


hunting <- rbind(N_bird, N_hunt)

SGN <- subset(hunting, grepl('SAGR', Species))

dodge = position_dodge(.5)
ggplot(SGN, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Sage-grouse Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

CHN <- subset(hunting, grepl('CHUK', Species))

dodge = position_dodge(.5)
ggplot(CHN, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Chukar Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

BLG <- subset(hunting, grepl('BLGR', Species))

dodge = position_dodge(.5)
ggplot(BLG, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Blue-grouse Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

HUP <- subset(hunting, grepl('HUPA', Species))

dodge = position_dodge(.5)
ggplot(HUP, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Hungarian Partridge Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

CAQ <- subset(hunting, grepl('CAQU', Species))

dodge = position_dodge(.5)
ggplot(CAQ, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and California Quail Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')

PHE <- subset(hunting, grepl('PHEA', Species))

dodge = position_dodge(.5)
ggplot(PHE, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region)) + #, formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Ring-necked Pheasant Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')


RAB <- subset(hunting, grepl('RABB', Species))

dodge = position_dodge(.5)
ggplot(RAB, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Region), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color =Region)) + #, formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) + labs(y = 'Trends in Hunter Numbers and Rabbit Harvested') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y') + theme_modern(legend.position = 'top')


lambda_bird <- MCMCsummary(mcmcList2, 'lambda')
lambda_hunt <- MCMCsummary(mcmcList2, 'lambda1')
lambda_survey <- MCMCsummary(mcmcList2, 'lambda2')


species_mod <- species[-c(4,8)]

lambda_bird$Species <- rep(species_mod, times = 41 * 2)
lambda_bird$Year   <- rep(1977:2017, each = 7 * 2)
lambda_bird$Region <- rep(c("West",'East'), each = 7, times = 41)
lambda_bird$Type <- rep('Bird')

ggplot(lambda_bird, aes(x = Year, y = `50%`, color = Region)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, fill = Region), shape  =21, color = 'black') +
  scale_fill_flat_d() +
  facet_wrap(~Species, ncol = 2, scales = 'free_y') + theme_bw() + theme_modern()

lambda_hunt$Species <- rep(species_mod, times = 41 * 2)
lambda_hunt$Year   <- rep(1977:2017, each = 7 * 2)
lambda_hunt$Region <- rep(c("West",'East'), each = 7, times = 41)
lambda_hunt$Type <- rep('Hunter')

ggplot(lambda_hunt, aes(x = Year, y = log(`50%`), color = Region)) +
  geom_pointrange(aes(x = Year, ymin = log(`2.5%`), ymax = log(`97.5%`), fill = Region), shape  =21, color = 'black') +
  scale_fill_flat_d() +
  facet_wrap(~Species, ncol = 3, scales = 'free_y') + theme_bw() + theme_modern()



lambdas <- rbind(lambda_bird, lambda_hunt)

SG <- subset(lambdas, grepl('SAGR', Species))

dodge = position_dodge(.5)
ggplot(SG, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()


CH <- subset(lambdas, grepl('CHUK', Species))
ggplot(CH, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','black')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()


CQ <- subset(lambdas, grepl('CAQU', Species))
ggplot(CQ, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','black')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()


HP <- subset(lambdas, grepl('HUPA', Species))
ggplot(HP, aes(x = Year, y = `50%`, color = Type)) +
  geom_pointrange(aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, color = Type), position = dodge, linetype = 'dashed' ) +
  geom_smooth( aes(x= Year, y = `50%`, color = Type), formula = y~x, method = 'lm') +
  geom_hline(yintercept = 1, color = 'red') +
  scale_color_manual(values = c('dodgerblue','black')) +
  facet_wrap(~Region, ncol = 1, scales = 'free_y') + theme_bw()