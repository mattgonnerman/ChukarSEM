### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Economic Predictors ###
  for(i in 1:4){
    sig.econ.pred[i] ~ dgamma(1,1)
    mu.econ.p[i] ~ dnorm(0, 1)
  }
  zeta.econ[1] <- 1
  for(i in 2:4){
    zeta.econ[i] ~ dnorm(0, 1)
  }
  for(t in 1:n.year){
    pred.econ.prime[t] ~ dnorm(0, 1) # Latent Predator Index
    gas[t] ~ dnorm(mu.econ.pred[1,t], sd = sig.econ.pred[1])
    une[t] ~ dnorm(mu.econ.pred[2,t], sd = sig.econ.pred[2])
    res[t] ~ dnorm(mu.econ.pred[3,t], sd = sig.econ.pred[3])
    pdi[t] ~ dnorm(mu.econ.pred[4,t], sd = sig.econ.pred[4])
    
    mu.econ.pred[1,t] <- mu.econ.p[1] + zeta.econ[1] * pred.econ.prime[t]
    mu.econ.pred[2,t] <- mu.econ.p[2] + zeta.econ[2] * pred.econ.prime[t]
    mu.econ.pred[3,t] <- mu.econ.p[3] + zeta.econ[3] * pred.econ.prime[t]
    mu.econ.pred[4,t] <- mu.econ.p[4] + zeta.econ[4] * pred.econ.prime[t]
  }
  
  ### Predator Predictor ###
  for(i in 1:4){
    sig.bbs.pred[i] ~ dgamma(1,1)
    mu.bbs.p[i] ~ dnorm(0, 1)
  }
  zeta.bbs[1] <- 1
  for(i in 2:4){
    zeta.bbs[i] ~ dnorm(0, 1)
  }
  for(t in 1:n.year){
    pred.bbs.prime[t] ~ dnorm(0, 1) # Latent Predator Index
    ravens[t] ~ dnorm(mu.bbs.pred[1,t], sd = sig.bbs.pred[1])
    rthawk[t] ~ dnorm(mu.bbs.pred[2,t], sd = sig.bbs.pred[2])
    nharr[t] ~ dnorm(mu.bbs.pred[3,t], sd = sig.bbs.pred[3])
    pfal[t] ~ dnorm(mu.bbs.pred[4,t], sd = sig.bbs.pred[4])
    
    mu.bbs.pred[1,t] <- mu.bbs.p[1] + zeta.bbs[1] * pred.bbs.prime[t]
    mu.bbs.pred[2,t] <- mu.bbs.p[2] + zeta.bbs[2] * pred.bbs.prime[t]
    mu.bbs.pred[3,t] <- mu.bbs.p[3] + zeta.bbs[3] * pred.bbs.prime[t]
    mu.bbs.pred[4,t] <- mu.bbs.p[4] + zeta.bbs[4] * pred.bbs.prime[t]
  }

  #Winter Severity
  beta.awssi ~ dnorm(0, sd = 5)
  alpha.awssi ~ dnorm(0, sd = 5)
  for (r in 1:n.region) {
    sig.awssi[r] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for (t in 1:n.year) {
      awssi[t,r] ~ dnorm(alpha.awssi + (beta.awssi * era.awssi[t]), sd = sig.awssi[r])
    }
  }
  
  # Drought Index
  for(r in 1:n.region){
    sig.drought[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(t in 1:n.year){
      pdsi[t,r] ~ dnorm(0, sd = sig.drought[r])
    } #t
  } #r
  
  #Bobcat Productivity
  sig.bob ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(t in 1:n.year){
    bobcat[t] ~ dnorm(mean = 0, sd = sig.bob)
  }
  
  
  
  ################################################################################
  ### Hunter Effort ###
  mu.econ ~ dnorm(0, sd = 100)
  sig.econ ~ T(dt(0, pow(2.5, -2), 1), 0, )
  for (s in 1:n.species) {
    beta.econ.hunt[s] ~ dnorm(mu.econ, sd = sig.econ)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      sig.H[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      alpha.hunt[s,r] ~ dnorm(mean = 0, sd = 7)
      for(k in 1:K){
        beta.spl.hunt[s,r,k] ~ dnorm(0, sd = sig.spl.hunt[s,r])
      } #k
      sig.spl.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(t in 1:n.year){ 
        
        latent.trend[s,t,r] <- inprod(beta.spl.hunt[s,r,1:K], basis[t,1:K]) #spline smoothing
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,t,r] <- alpha.hunt[s,r] + #intercept
          beta.econ.hunt[s] * pred.econ.prime[t] + #SEM economic indicator
          latent.trend[s,t,r]
        
        H[s,t,r] <- exp(hunt.eps[s,t,r]) #Log Link
        
        n.hunt[s,t,r] ~ dnorm(H[s,t,r], sd = sig.H[s,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.hunt[s,t,r] <- H[s,t+1,r]/H[s,t,r]
        log.r.hunt[s,t,r] <- log(lambda.hunt[s,t,r])
      } #t
    } #s
    
    for(t in 1:n.year){
      hunt.eps[1:n.species,t,r] ~ dmnorm(mu.hunt[1:n.species,t,r], cov =  Sigma.hunt[1:n.species,1:n.species,r] )
    } #t
    
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
  mu.bobcat ~ dnorm(0, 0.01)
  sig.bobcat ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.bbs ~ dnorm(0, 0.01)
  sig.bbs ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.pdsi ~ dnorm(0, 0.01)
  sig.pdsi ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
      beta.wintsev.harv[s] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
      beta.bbs.harv[s] ~ dnorm(mu.bbs, sd  = sig.bbs)
      beta.pdsi.harv[s] ~ dnorm(mu.pdsi, sd  = sig.pdsi)
      beta.bobcat.harv[s] ~ dnorm(mu.bobcat, sd  = sig.bobcat)
  } #s
  
  for(r in 1:n.region){
    mu.hunter.harv[r] ~ dnorm(0, 0.01)
    sig.hunter.harv[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for(s in 1:n.species){
      beta.hunter.harv[s,r] ~ dnorm(mu.hunter.harv[r], sd = sig.hunter.harv[r])
      sig.N[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      alpha.harv[s,r] ~ dnorm(mean = 0, sd = 7)
      
      # Process Model
      for(t in 1:n.year){
        #Unlinked estimate of Hunter Numbers
        mu.harv[s,t,r] <- alpha.harv[s,r] + # intercepts
          beta.hunter.harv[s,r] * latent.trend[s,t,r] + # Latent Hunter Trend
          beta.wintsev.harv[s] * awssi[t,r]  + # Previous winter severity (Affecting Survival)
          # beta.pdsi.harv[s] * pdsi[t,r] + # Same year, spring/summer drought (Affecting Survival/Reproduction)
          beta.bbs.harv[s] * pred.bbs.prime[t] #+ # Latent predator index (Affecting Reproduction)
          # beta.bobcat.harv[s] * bobcat[t]
        
        N[s,t,r] <- exp(harv.eps[s,t,r]) #Log Link
        
        n.harv[s,t,r] ~ dnorm(N[s,t,r], sd = sig.N[s,r]) #Number of hunters follows Poisson
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
    mu.chuk[r] ~ dlogis(0,1)
  }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.chuk){
      log.r.chuk[p,t-1] <-  mod.chuk[reg.chuk[p]] * log.r.harv[3, t+13, reg.chuk[p]] 
      
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
        BPH[s,t,r] <-  N[s,t,r]/H[s,t,r]
        BPH2[s,t,r]<- exp(mu.hunt[s,t,r]-mu.harv[s,t,r])
      }
    }
  }
   
  # for(s in 1:n.species){
  #   sig.bph[s] ~ T(dt(0, pow(2.5,-2), 1),0,)
  #   for(t in 1:3){
  #     bph.survey[s,t] ~ dnorm(mean = mean(BPH[s,42+t,1:2]), sd = sig.bph[s])
  #   }
  # }
})