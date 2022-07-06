### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Predictors ###
  #Ravens (Highly correlated with Prairie Falcon and RTHawk)
  alpha.rav ~ dnorm(0, 0.001) #intercept
  beta.t.rav ~ dnorm(0, 0.01) #year
  sig.rav~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1.rav ~ dunif(-1,1) #Autoregressive parameter
  
  rav.trend[1] <- alpha.rav + beta.t.rav * 1
  mu.rav[1] <- rav.trend[1]
  for(t in 2:n.year){
    rav.trend[t] <- alpha.rav + beta.t.rav * t
    mu.rav[t] <- rav.trend[t] + ar1.rav * (raven[t-1] - rav.trend[t-1])
  } #t
  
  for(t in 1:n.year){
    raven[t] ~ dnorm(mu.rav[t], sd = sig.rav)
  }
  
  #Unemployment (Jobs)
  sig.une ~ T(dt(0, pow(2.5, -2), 1), 0, )
  for (t in 1:(n.year)) {
    une[t] ~ dnorm(0, sd = sig.une)
  }
  
  #Winter Severity
  beta.awssi ~ dnorm(0, sd = 100)
  alpha.awssi ~ dnorm(0, sd = 100)
  for (r in 1:n.region) {
    sig.awssi[r] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for (t in 1:(n.year+1)) {
      awssi[r, t] ~ dnorm(alpha.awssi + beta.awssi * (era.awssi[t] - 1), sd = sig.awssi[r])
    }
  }
  
  ################################################################################
  ### Hunter Effort ###
  mu.wintsev.hunt ~ dnorm(0, sd = 100)
  sig.wintsev.hunt ~ T(dt(0, pow(2.5, -2), 1), 0, )
  mu.jobs ~ dnorm(0, sd = 100)
  sig.jobs ~ T(dt(0, pow(2.5, -2), 1), 0, )
  
  for (s in 1:n.species) {
    beta.wintsev.hunt[s] ~ dnorm(mu.wintsev.hunt, sd = sig.wintsev.hunt)
    beta.jobs[s] ~ dnorm(mu.jobs, sd = sig.jobs)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.hunt[s,r] ~ dnorm(5, sd = 12) #sd = 1)
      for(k in 1:K){
        beta.spl.hunt[s,r,k] ~ dnorm(0, sd = sig.spl.hunt[s,r])
      } #k
      sig.spl.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(t in 1:n.year){ 
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,t,r] <- alpha.hunt[s,r] + #intercept
          beta.wintsev.hunt[s] * awssi[r, t + 1] + 
          beta.jobs[s] * une[t] + 
          inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #spline smoothing
        
        H[s,t,r] <- exp(hunt.eps[s,t,r]) #Log Link
        
        n.hunt[s,t,r] ~ dpois(H[s,t,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.hunt[s,t,r] <- H[s,t+1,r]/H[s,t,r]
        log.r.hunt[s,t,r] <- log(lambda.hunt[s,t,r])
      }
    } #s
    
    for(t in 1:n.year){
      hunt.eps[1:n.species,t,r] ~ dmnorm(mu.hunt[1:n.species,t,r], cov =  Sigma.hunt[1:n.species,1:n.species,r] )
    }
    
    # Correlation Matrices
    Q.hunt[1:n.species,1:n.species,r] ~ dinvwish(S = I.hunt[1:n.species,1:n.species,r], df = n.species + 1)
    
    for(s in 1:n.species){
      sig.hunt[s,r] ~ dgamma(1,1)
      Delta.hunt[s,s,r] <- pow(Q.hunt[s,s,r], -0.5)
      Lambda.hunt[s,s,r] <- sig.hunt[s,r]
    } #s
    
    for (s1 in 2:n.species){
      for (s2 in 1:(s1-1)){
        Lambda.hunt[s1,s2,r] <- 0
        Delta.hunt[s1,s2,r] <- 0
      } #s2
    } #s1
    
    Sigma.hunt[1:n.species,1:n.species,r] <- Lambda.hunt[1:n.species,1:n.species,r] %*% P.hunt[1:n.species,1:n.species,r] %*% Lambda.hunt[1:n.species,1:n.species,r]  
    P.hunt[1:n.species,1:n.species,r] <- Delta.hunt[1:n.species,1:n.species,r] %*% Q.hunt[1:n.species,1:n.species,r] %*% Delta.hunt[1:n.species,1:n.species,r]
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho.hunt[s1,s2,r] <- Sigma.hunt[s1,s2,r]/sqrt(Sigma.hunt[s1,s1,r] * Sigma.hunt[s2,s2,r])   
      } #s2
    } #s1
  } #r
  
  
  ################################################################################
  ### Total Harvest ###
  mu.wintsev.harv ~ dnorm(0, 0.01)
  sig.wintsev.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.raven ~ dnorm(0, 0.01)
  sig.raven ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
    beta.wintsev.harv[s] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
    beta.raven.harv[s] ~ dnorm(mu.raven, sd  = sig.raven)
  } #s
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.harv[s,r] ~ dnorm(0, sd = 12)
      for(k in 1:K){
        beta.spl.harv[s,r,k] ~ dnorm(0, sd = sig.spl.harv[s,r])
      } #k
      sig.spl.harv[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      for(t in 1:n.year){
        #Unlinked estimate of Hunter Numbers
        mu.harv[s,t,r] <- alpha.harv[s,r] +#regression formula
          beta.wintsev.harv[s] * awssi[r,t] + #Previous winter severity (Affecting Survival)
          beta.raven.harv[s] * raven[t] + #Previous Year BBS index (Affecting Reproduction)
          inprod(beta.spl.harv[s,r,1:K], Z.harv[t,1:K,s,r]) #spline smoothing
        
        N[s,t,r] <- exp(harv.eps[s,t,r]) #Log Link
        
        n.harv[s,t,r] ~ dpois(N[s,t,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.harv[s,t,r] <- N[s,t+1,r]/N[s,t,r]
        log.r.harv[s,t,r] <- log(lambda.harv[s,t,r])
      }
    } #s
    
    for(t in 1:n.year){
      harv.eps[1:n.species,t,r] ~ dmnorm(mu.harv[1:n.species,t,r], cov =  Sigma.harv[1:n.species,1:n.species,r] )
    }
    
    ### Correlation Matrices
    Q.harv[1:n.species,1:n.species,r] ~ dinvwish(S = I.harv[1:n.species,1:n.species,r], df = n.species + 1)
    
    for(s in 1:n.species){
      sig.harv[s,r] ~ dgamma(1,1)
      Delta.harv[s,s,r] <- pow(Q.harv[s,s,r], -0.5)
      Lambda.harv[s,s,r] <- sig.harv[s,r]
    } #s
    
    for (s1 in 2:n.species){
      for (s2 in 1:(s1-1)){
        Lambda.harv[s1,s2,r] <- 0
        Delta.harv[s1,s2,r] <- 0
      } #s2
    } #s1
    
    Sigma.harv[1:n.species,1:n.species,r] <- Lambda.harv[1:n.species,1:n.species,r] %*% P.harv[1:n.species,1:n.species,r] %*% Lambda.harv[1:n.species,1:n.species,r]  
    P.harv[1:n.species,1:n.species,r] <- Delta.harv[1:n.species,1:n.species,r] %*% Q.harv[1:n.species,1:n.species,r] %*% Delta.harv[1:n.species,1:n.species,r]  
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho.harv[s1,s2,r] <- Sigma.harv[s1,s2,r]/sqrt(Sigma.harv[s1,s1,r] * Sigma.harv[s2,s2,r])   
      } #s2
    } #s1
  } #r
  
  ################################################################################
  ### Chukar Site Abundance ###
  for(r in 1:n.region){
    theta.chuk[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    mod.chuk[r] ~ dlogis(0,1)
  }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.chuk){
      log.r.chuk[p,t-1] <-  mod.chuk[reg.chuk[p]] * log.r.harv[3, t+13, reg.chuk[p]] #log.r.harv[t=13] is 1990-1991
      
      C.chuk[p,t] <- exp(log.r.chuk[p,t-1]) * C.chuk[p,t-1] #Equivalent of Poisson lambda
      
      rate.chuk[p,t-1] <- theta.chuk[reg.chuk[p]]/(theta.chuk[reg.chuk[p]] + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(prob = rate.chuk[p,t-1], size = theta.chuk[reg.chuk[p]]) #obs. # of chukars follow neg-bin
    } #t
  } #p
  
  ################################################################################
  ### Birds per Hunter
  for(t in 1:n.year){
    for(s in 1:n.species){
      for(r in 1:n.region){
        BPH[s,t,r] <- N[s,t,r]/H[s,t,r]
      }
    }
  }
})