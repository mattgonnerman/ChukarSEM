require(dplyr)

remove('n.bias', 'h.bias', 'bph.bias')
for(i in 2012:2017){
  full.df <- read.csv(paste("./Output/Model Validation/CHUK - ",i, "/EstCheck - ChukOnly - Estimates_N_H_BPH.csv", sep = "")) %>%
    filter(Year %in% (i):(i+2)) %>%
    filter(!is.na(Obs.BPHS))
  
  n.df <- full.df %>% select(County, Year, Est = Est.N, LCL = Est.N.LCL, UCL = Est.N.UCL, Obs = Obs.N) %>%
    filter(!is.na(Obs))%>%
    mutate(Est = Est*1000,
           LCL = LCL*1000,
           UCL = UCL*1000) %>%
    mutate(InOut = ifelse(Obs >= LCL & Obs <= UCL, 1, 0),
           RelBias = (Obs - Est)/Obs) %>%
    select(County, Year, InOut, RelBias)
  
  h.df <- full.df %>% select(County, Year, Est = Est.H, LCL = Est.H.LCL, UCL = Est.H.UCL, Obs = Obs.H) %>%
    filter(!is.na(Obs))%>%
    mutate(Est = Est*1000,
           LCL = LCL*1000,
           UCL = UCL*1000) %>%
    mutate(InOut = ifelse(Obs >= LCL & Obs <= UCL, 1, 0),
           RelBias = (Obs - Est)/Obs) %>%
    select(County, Year, InOut, RelBias)
  
  bph.df <- full.df %>% select(County, Year, Est = Est.BPH, LCL = Est.BPH.LCL, UCL = Est.BPH.UCL, Obs = Obs.BPHS) %>%
    filter(!is.na(Obs))%>%
    mutate(Est = Est,
           LCL = LCL,
           UCL = UCL) %>%
    mutate(InOut = ifelse(Obs >= LCL & Obs <= UCL, 1, 0),
           RelBias = (Obs - Est)/Obs) %>%
    select(County, Year, InOut, RelBias)
  
  if(!exists('n.bias')){
    n.bias <- n.df
    h.bias <- h.df
    bph.bias <- bph.df
  }else{
    n.bias <- rbind(n.bias, n.df)
    h.bias <- rbind(h.bias, h.df)
    bph.bias <- rbind(bph.bias, bph.df)
  }
}

length(n.bias[n.bias$InOut == 1,]$InOut)
length(n.bias$InOut)
summary(n.bias$RelBias)

length(h.bias[h.bias$InOut == 1,]$InOut)
length(h.bias$InOut)
summary(h.bias$RelBias)

length(bph.bias[bph.bias$InOut == 1,]$InOut)
length(bph.bias$InOut)
summary(bph.bias$RelBias)

