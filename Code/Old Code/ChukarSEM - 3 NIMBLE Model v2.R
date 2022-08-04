code <- nimbleCode( {
  for(r in 1:n.region){
    sig.drought[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.drought2[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.jobs[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    
    for(s in 1:n.species){
      lbo1[s,r] ~ dnorm(5, sd = 1) #Intercept for Number of Hunters
      lbo[s,r] ~ dnorm(0, sd = 1) #Intercept for Number of Birds Harvested
      beta.drought[s,r] ~ dnorm(0, sd = sig.drought[r]) # Effect of drought on 
      beta.drought2[s,r] ~ dnorm(0, sd = sig.drought2[r]) 
      beta.jobs[s,r] ~ dnorm(0, sd = sig.jobs[r]) 
      
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
        
        #Unlinked estimate of Hunter Numebrs
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
      
      #Track Change in Harvest
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Number Harvested, Year 1
      for(t in 2:(n.year)){
        #Spline trend
        pred[s,t-1,r] <- inprod(beta.pressure[s,r,1:K], Z[t-1,1:K,s,r])
        
        #Regression for change in species-specific total harvest
        mu2[s,t-1,r] <- lbo[s,r] + 
          inprod(beta.pressure[s,r,1:K], Z[t-1,1:K,s,r]) + 
          beta.drought[s,r] * pdsi[t-1,r]
        
        lambda[s,t-1,r] <- exp(log.r[s,t-1,r]) #link function for growth rate
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
      C[p,1] ~ dpois(x[p,1]) #Equivalent of Poisson lambda
      sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
      mod[p] ~ dlogis(0,1)
      
      for(t in 1:(n.yr-1)){
        chuk.eps[p,t]  ~ dnorm(mod[p] * log.r[3, t+14, 1], sd = sigma.chuk[p]) #Coeffient for annual change
        log.r2[p,t] <- lbo2[p] + chuk.eps[p,t] #Unlinked change in abundance
        lambda2[p,t] <- exp(log.r2[p,t]) #Change in abundance
        C[p,t+1] <- lambda2[p,t] * C[p,t] #Equivalent of Poisson lambda
        rate2[p,t+1] <- theta2[p]/(theta2[p] + C[p,t+1]) #NB success parameter
        n.chuk[p,t+1] ~ dnegbin(rate2[p,t+1], theta2[p]) #obs. # of chukars follow neg-bin  
      } #t
  } #p  
    
  ### Correlation Matrices
  # from Barnard et al. 2000. Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage
  # Why we used it rather than standard https://arxiv.org/pdf/1408.4050.pdf
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
    
    #Correlation Coefficients
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho[s1,s2,r] <- Sigma[s1,s2,r]/sqrt(Sigma[s1,s1,r] * Sigma[s2,s2,r]) #Number of Hunters
        rho2[s1,s2,r] <- Sigma2[s1,s2,r]/sqrt(Sigma2[s1,s1,r] * Sigma2[s2,s2,r]) #Total Harvest
      } #s2
    } #s1
  } #r

})
