code <- nimbleCode( {
  #Unemployment
  sig.une ~ dunif(0,5)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une) 
  } #t
  
  for(r in 1:n.region){
    #Drought Conditions
    sig.wpdsi[r] ~ dunif(0,5)
    sig.pdsi[r] ~ dunif(0,5)
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
      pdsi[t,r] ~ dnorm(0, sd = sig.pdsi[r])
    } #t       
    
    for(s in 1:n.species){
      lbo1[s,r] ~ dnorm(5, sd = 1)
      lbo[s,r] ~ dnorm(0, sd = 1)
      beta.drought[s,r] ~ dnorm(mu.drought[r], sd = sig.drought[r])
      beta.drought2[s,r] ~ dnorm(mu.drought2[r], sd = sig.drought2[r]) 
      beta.jobs[s,r] ~ dnorm(mu.jobs[r], sd = sig.jobs[r]) 
      
      #For Splines
      sig.pressure[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      sig.trend[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      for(k in 1:K){
        beta.pressure[s,r,k] ~ dnorm(0, sd = sig.pressure[s,r]) 
        beta.trend[s,r,k] ~ dnorm(0, sd = sig.trend[s,r])
      } #k
      
      for(t in 1:(n.year)){
        #Spline Trends
        pred1[s,r,t] <- inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r])
        
        #Unlinked estimate of Hunter Numbers
        mu[s,r,t] <- lbo1[s,r] + 
          inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r]) +
          beta.drought2[s,r] * wpdsi[t,r] +
          beta.jobs[s,r] * une[t]
        hunt.eps[1:n.species,r,t] ~ dmnorm( mu[1:n.species,r,t], cov =  Sigma[1:n.species,1:n.species,r] )
        log(H[s,t,r]) <- hunt.eps[s,r,t] #Log Link
        n.hunter[s,t,r] ~  dpois(H[s,t,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 2:(n.year)){
        lambda1[s,t-1,r] <- H[s,t,r]/H[s,t-1,r] #Change in hunters
      } #t
      
      #Track Change in Total Harvest
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Total harvest, Year 1
      for(t in 2:(n.year)){
        #Spline trend
        pred[s,t-1,r] <- inprod(beta.pressure[s,r,1:K], Z[t-1,1:K,s,r])
        
        #Regression for change in total harvest
        mu2[s,t-1,r] <- lbo[s,r] + 
          inprod(beta.pressure[s,r,1:K], Z[t-1,1:K,s,r]) + 
          beta.drought[s,r] * pdsi[t-1,r]
        
        lambda[s,t-1,r] <- exp(log.r[s,t-1,r]) #link function
        N[s,t,r] <- lambda[s,t-1,r] * N[s,t-1,r] #number available = change since last year
        n.harv[s,t,r] ~  dpois(N[s,t,r]) #Number harvested follows Poisson
        
        bph[s,t,r] <-  N[s,t,r]/H[s,t,r] #birds per hunter
      } #t
      bph[s,1,r] <-  N[s,1,r]/H[s,1,r] #birds per hunter, Year 1
    } #s 
    
    for(t in 2:(n.year)){
      #Change in total harvest
      log.r[1:n.species,t-1,r]  ~ dmnorm( mu2[1:n.species,h-1,r], cov =  Sigma2[1:n.species,1:n.species,r])
    } #t
    
  } #r

  
  ### Chukar Site Abundance
    for(p in 1:n.site){
      lbo2[p] ~ dlogis(0,1) #Intercept
      theta2[p] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB overdispersion parameter
      C[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
      sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
      mod[p] ~ dlogis(0,1)
      
      for(t in 2:n.yr){
        chuk.eps[p,t-1]  ~ dnorm(mod[p] * log.r[3, t+13, 1], sd = sigma.chuk[p]) #Coeffient for annual change
        log.r2[p,t-1] <- lbo2[p] + chuk.eps[p,t-1] #Unlinked change in abundance
        lambda2[p,t-1] <- exp(log.r2[p,t-1]) #Change in abundance
        C[p,t] <- lambda2[p,t-1] * C[p,t-1] #Equivalent of Poisson lambda
        rate2[p,t] <- theta2[p]/(theta2[p] + C[p,t]) #NB success parameter
        n.chuk[p,t] ~ dnegbin(rate2[p,t], theta2[p]) #obs. # of chukars follow neg-bin
      } #t
    } #p  
    
  ### Correlation Matrices
  for(r in 1:n.region){
    Q[1:n.species,1:n.species,r] ~ dinvwish(S = I[1:n.species,1:n.species,r], df = n.species + 1)
    Q2[1:n.species,1:n.species,r] ~ dinvwish(S = I2[1:n.species,1:n.species,r], df = n.species + 1)
    
    for(s in 1:n.species){
      sig[s,r] ~ dgamma(1,1)
      Delta[s,s,r] <- pow(Q[s,s,r], -0.5)
      Lambda[s,s,r] <- sig[s,r]
      
      sig2[s,r] ~ dgamma(1,1)
      Delta2[s,s,r] <- pow(Q2[s,s,r], -0.5)
      Lambda2[s,s,r] <- sig2[s,r]
    } #s
    
    for (s1 in 2:n.species){
      for (s2 in 1:(s1-1)){
        Lambda[s1,s2,r] <- 0
        Delta[s1,s2,r] <- 0
        
        Lambda2[s1,s2,r] <- 0
        Delta2[s1,s2,r] <- 0
      } #s2
    } #s1
    
    Sigma[1:n.species,1:n.species,r] <- Lambda[1:n.species,1:n.species,r] %*% P[1:n.species,1:n.species,r] %*% Lambda[1:n.species,1:n.species,r]  
    P[1:n.species,1:n.species,r] <- Delta[1:n.species,1:n.species,r] %*% Q[1:n.species,1:n.species,r] %*% Delta[1:n.species,1:n.species,r]
    
    Sigma2[1:n.species,1:n.species,r] <- Lambda2[1:n.species,1:n.species,r] %*% P2[1:n.species,1:n.species,r] %*% Lambda2[1:n.species,1:n.species,r]  
    P2[1:n.species,1:n.species,r] <- Delta2[1:n.species,1:n.species,r] %*% Q2[1:n.species,1:n.species,r] %*% Delta2[1:n.species,1:n.species,r]  
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho[s1,s2,r] <- Sigma[s1,s2,r]/sqrt(Sigma[s1,s1,r] * Sigma[s2,s2,r])   
        rho2[s1,s2,r] <- Sigma2[s1,s2,r]/sqrt(Sigma2[s1,s1,r] * Sigma2[s2,s2,r])   
      } #s2
    } #s1
  } #r

})
