code <- nimbleCode( {
  ### Economic Trends
  #Personal Disposable Income
  b0.pdi ~ dnorm(0, 0.001) #intercept
  bt.pdi ~ dnorm(0, 0.01) #year
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1 ~ dunif(-1,1) #Change in PDI over time?
  
  pdi.trend[1] <- b0.pdi + bt.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    pdi.trend[t] <- b0.pdi + bt.pdi * t
    mu.pdi[t] <- pdi.trend[t] + ar1 * (PDI[t-1] - pdi.trend[t-1])
  } #t
  
  #Gas Prices
  b0.gas ~ dnorm(1.5, 1)
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  bt.gas[1] ~ dnorm(0, 0.01)
  bt.gas[2] ~ dnorm(0, 0.01)
  
  #Relative Cost of Gas
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(b0.gas + bt.gas[era[t]]*t, sd = sig.gas),0,)
    REL.COST[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  } #t
  
  #Unemployment Rate
  sig.une ~ dunif(0,5)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
  } #t
  
  ### Trends in Harvest
  for(r in 1:n.region){
    sig.wpdsi[r] ~ dunif(0,5)
    
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
      hunt.eps[1:n.species,r,t] ~ dmnorm( mu[1:n.species,r,t], cov =  Sigma[1:n.species,1:n.species,r] )
    } #t
  } #r
  
  for (s in 1:n.species){
    for(r in 1:n.region){
      beta.drought2[s,r] ~ dnorm(mu.drought2[r], sd = sig.drought2[r]) 
      beta.jobs[s,r] ~ dnorm(mu.jobs[r], sd = sig.jobs[r]) 
      beta.general[s,r] ~ dnorm(mu.gen[r], sd  = sig.gen[r])
      beta.income[s,r] ~ dnorm(mu.incom[r], sd  = sig.incom[r])
      
      lbo1[s,r] ~ dnorm(5, sd = 1)
      
      log(H[s,1,r]) <- hunt.eps[s,r,1]
      z[s,1,r] ~  dpois(H[s,1,r]) #number of hunters observed
      for(t in 2:(n.year)){
        log(H[s,t,r]) <- hunt.eps[s,r,t]
        z[s,t,r] ~  dpois(H[s,t,r])
      } #t
    } #r
  } #s
  
  for(r in 1:n.region){
    mu.incom[r] ~ dnorm(0, 0.01)
    sig.incom[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.gen[r] ~ dnorm(0, 0.01)
    sig.gen[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.drought2[r] ~ dnorm(0, 0.01)
    sig.drought2[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.jobs[r] ~ dnorm(0, 0.01)
    sig.jobs[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    
    for (s in 1:n.species){
      sig.trend[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(k in 1:K){
        beta.trend[s,r,k] ~ dnorm(0, sd = sig.trend[s,r])
      } #k
    } #s
    
    for(t in 2:(n.year)){
      for (s in 1:n.species){
        lambda1[s,t-1,r] <- H[s,t,r]/H[s,t-1,r]
      } #s
    } #t
  } #r
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      sig[s,r] ~ dgamma(1,1)
      Delta[s,s,r] <- pow(Q[m,m,r], -0.5)
      Lambda[s,s,r] <- sig[m,r]
    } #s
    
    for (s1 in 2:n.species){
      for (s2 in 1:(i-1)){
        Lambda[s1,s2,r] <- 0
        Delta[s1,s2,r] <- 0
      } #j
    } #i
    
    #
    Q[1:n.species,1:n.species,r] ~ dinvwish(S = I[1:n.species,1:n.species,r], df = n.species + 1)
    P[1:n.species,1:n.species,r] <- Delta[1:n.species,1:n.species,r] %*% Q[1:n.species,1:n.species,r] %*% Delta[1:n.species,1:n.species,r]
    Sigma[1:n.species,1:n.species,r] <- Lambda[1:n.species,1:n.species,r] %*% P[1:n.species,1:n.species,r] %*% Lambda[1:n.species,1:n.species,r]  
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho[s1,s2,r] <- Sigma[s1,s2,r]/sqrt(Sigma[s1,s1,r] * Sigma[s2,s2,r])   
      } #s2
    } #s1
    
    for(s in 1:n.species){
      for(t in 1:n.year){
        pred1[s,r,t] <- inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r])
        mu[s,r,t] <- lbo1[s,r] + #
          inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r]) + #Spline smoothing terms
          beta.general[s,r] * res[t] + #Hunting Licences Sold
          beta.income[s,r] * REL.COST[t] + #PDI/Gas Price
          beta.drought2[s,r] * wpdsi[t,r] + #Drought
          beta.jobs[s,r] * une[t] #Unemployment
      } #t
    } #s
  } #r
  
})