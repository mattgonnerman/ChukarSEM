library(coda)
library(MCMCvis)

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

test1 <- MCMCsummary(mcmcList2, 'pred1')
test1.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)],  region = c('east','west'), year = 1976:2021 ), test1)

fig1 <- ggplot(data = test1.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),
                  shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) +
  labs(y = 'Trends in Hunter Numbers', x = 'Year') +
  facet_grid(region ~ species)
  # scale_color_flat::scale_fill_flat_d() +
  # theme_modern()
fig1

test3 <- MCMCsummary(mcmcList2, 'pred2')
test3.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)],  region = c('east','west'), year = 1976:2021 ), test3)

fig3 <- ggplot(data = test3.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_hline(yintercept = 0) + labs(y = 'Trends in Hunter Numbers', x = 'Year') + facet_grid(region ~ species) +  scale_fill_flat_d() + theme_modern()
fig3



hunter_prime <- abind(hunters, array(NA, dim = c(7,4,2)), along = 2)
test2    <- MCMCsummary(mcmcList2, 'H')
test2.df <- cbind.data.frame(expand.grid(species = species[-c(4,8)], year = 1976:2021, region = c('east','west') ), test2, hunters = c(hunter_prime))

fig2 <- ggplot(data = test2.df, aes(x = year, y = mean, group = region)) +
  geom_pointrange2(aes(x =year, y = mean, ymin = `2.5%`,ymax = `97.5%`, fill = species),shape  = 21, color = 'black', size = .75, position = position_dodge(1)) +
  geom_point(aes(x =year, y = hunters )) +
  geom_hline(yintercept = 0) + labs(y = 'Hunter Numbers', x = 'Year') + facet_wrap(region ~ species, ncol  = 7, scales = 'free') +  scale_fill_flat_d() + theme_modern()
fig2