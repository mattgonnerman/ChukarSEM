code <- nimbleCode( {
  for(r in 1:n.region){
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
      }
    }  
    
    Sigma[1:n.species,1:n.species,r] <- Lambda[1:n.species,1:n.species,r] %*% P[1:n.species,1:n.species,r] %*% Lambda[1:n.species,1:n.species,r]  
    Q[1:n.species,1:n.species,r] ~ dinvwish(S = I[1:n.species,1:n.species,r], df = n.species + 1)
    P[1:n.species,1:n.species,r] <- Delta[1:n.species,1:n.species,r] %*% Q[1:n.species,1:n.species,r] %*% Delta[1:n.species,1:n.species,r]
    
    Sigma2[1:n.species,1:n.species,r] <- Lambda2[1:n.species,1:n.species,r] %*% P2[1:n.species,1:n.species,r] %*% Lambda2[1:n.species,1:n.species,r]  
    Q2[1:n.species,1:n.species,r] ~ dinvwish(S = I2[1:n.species,1:n.species,r], df = n.species + 1)
    P2[1:n.species,1:n.species,r] <- Delta2[1:n.species,1:n.species,r] %*% Q2[1:n.species,1:n.species,r] %*% Delta2[1:n.species,1:n.species,r]  
    
    for (s1 in 1:n.species){
      for (s2 in 1:n.species){
        rho[s1,s2,r] <- Sigma[s1,s2,r]/sqrt(Sigma[s1,s1,r] * Sigma[s2,s2,r])   
        rho2[s1,s2,r] <- Sigma2[s1,s2,r]/sqrt(Sigma2[s1,s1,r] * Sigma2[s2,s2,r])   
      } #s2
    } #s1
    
    #Variation in Drought Conditions
    sig.wpdsi[r] ~ dunif(0,5)
    sig.pdsi[r] ~ dunif(0,5)
    #Year specific drought conditions
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
      pdsi[t,r] ~ dnorm(0, sd = sig.pdsi[r])
    }        
    
    for(s in 1:n.species){
      for(t in 1:n.year){
        pred1[s,r,t] <- inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r])
        mu[s,r,t] <- lbo1[s,r] + inprod(beta.trend[s,r,1:K], ZZ[t,1:K,s,r]) + beta.drought2[s,r] * wpdsi[t,r] + beta.jobs[s,r] * une[t]
      } #t
    } #s
  } #r
  
  sig.une~ dunif(0,5)
  
  #
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
    for(r in 1:n.region){
      hunt.eps[1:n.species,r,t] ~ dmnorm( mu[1:n.species,r,t], cov =  Sigma[1:n.species,1:n.species,r] )
    } #r
  } #t
  
  #
  for (s in 1:n.species){
    for(r in 1:n.region){
      beta.drought2[s,r] ~ dnorm(mu.drought2[r], sd = sig.drought2[r]) 
      beta.jobs[s,r] ~ dnorm(mu.jobs[r], sd = sig.jobs[r]) 
      lbo1[s,r] ~ dnorm(5, sd = 1)
      lbo[s,r] ~ dnorm(0, sd = 1)
      log(H[s,1,r]) <- hunt.eps[s,r,1]
      z[s,1,r] ~  dpois(H[s,1,r])
      
      for(t in 2:(n.year)){
        log(H[s,t,r]) <- hunt.eps[s,r,t]
        z[s,t,r] ~  dpois(H[s,t,r])
      }# t
    }
  
    
    mu.drought[r] ~ dnorm(0, 0.01)
    sig.drought[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.drought2[r] ~ dnorm(0, 0.01)
    sig.drought2[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.jobs[r] ~ dnorm(0, 0.01)
    sig.jobs[r] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for (i in 1:n.species){
      beta.drought[i,r] ~ dnorm(mu.drought[r], sd = sig.drought[r])
      sig.pressure[i,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      sig.trend[i,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      for(k in 1:K){
        beta.pressure[i,r,k] ~ dnorm(0, sd = sig.pressure[i,r]) 
        beta.trend[i,r,k] ~ dnorm(0, sd = sig.trend[i,r])
      }
    } 
    for (i in 1:n.species){
      N[i,1,r] ~ dpois(y[i,1,r])
      bph[i,1,r] <-  N[i,1,r]/H[i,1,r]
      for(h in 2:(n.year)){
        mu2[i,h-1,r] <- lbo[i,r] + inprod(beta.pressure[i,r,1:K], Z[h-1,1:K,i,r]) + beta.drought[i,r] * pdsi[h-1,r]
      }
    }
    for(h in 2:(n.year)){
      log.r[1:n.species,h-1,r]  ~ dmnorm( mu2[1:n.species,h-1,r], cov =  Sigma2[1:n.species,1:n.species,r])
    }
    for(h in 2:(n.year)){
      for (i in 1:n.species){
        pred[i,h-1,r] <- inprod(beta.pressure[i,r,1:K], Z[h-1,1:K,i,r]) 
        lambda[i,h-1,r] <- exp(log.r[i,h-1,r])
        lambda1[i,h-1,r] <- H[i,h,r]/H[i,h-1,r]
        N[i,h,r] <- lambda[i,h-1,r] * N[i,h-1,r]
        y[i,h,r] ~  dpois(N[i,h,r]) 
        bph[i,h,r] <-  N[i,h,r]/H[i,h,r]
      }
    }
  }
  
  for(p in 1:n.site){
    sigma.chuk[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
    lbo2[p] ~ dlogis(0,1)
    theta2[p] ~ T(dt(0, pow(2.5,-2), 1),0,)
    C[p,1] ~ dpois(x[p,1])
    mod[p] ~ dlogis(0,1)
    for(h in 2:n.yr){
      log.r2[p,h-1] <- lbo2[p] + chuk.eps[p,h-1]
      lambda2[p,h-1] <- exp(log.r2[p,h-1])
      C[p,h] <- lambda2[p,h-1] * C[p,h-1]
      rate2[p,h] <- theta2[p]/(theta2[p] + C[p,h])
      x[p,h] ~   dnegbin(rate2[p,h], theta2[p])  
    }
    for(h in 1:(n.yr-1)){
      chuk.eps[p,h]  ~ dnorm(mod[p] * log.r[3,h+14,1], sd = sigma.chuk[p])
    }
  }
  
})