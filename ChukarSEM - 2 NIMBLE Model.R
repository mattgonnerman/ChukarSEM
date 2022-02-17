code <- nimbleCode( {
  
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(b0.gas + bt.gas[era[t]]*t, sd = sig.gas),0,)
    REL.COST[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  }
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  b0.pdi ~ dnorm(0, 0.001)
  bt.pdi ~ dnorm(0, 0.01)
  b0.gas ~ dnorm(1.5, 1)
  bt.gas[1] ~ dnorm(0, 0.01)
  bt.gas[2] ~ dnorm(0, 0.01)
  ar1 ~ dunif(-1,1)
  pdi.trend[1] <- b0.pdi + bt.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    mu.pdi[t] <- pdi.trend[t] + ar1 * ( PDI[t-1] - pdi.trend[t-1] )
    pdi.trend[t] <- b0.pdi + bt.pdi * t
  }
  
  for(l in 1:n.region){
    for(m in 1:n.species){
      sig[m,l] ~ dgamma(1,1)
      Delta[m,m,l] <- pow(Q[m,m,l], -0.5)
      Lambda[m,m,l] <- sig[m,l]
    }
    for (i in 2:n.species){
      for (j in 1:(i-1)){
        Lambda[i,j,l] <- 0
        Delta[i,j,l] <- 0
      }
    }  
    Sigma[1:n.species,1:n.species,l] <- Lambda[1:n.species,1:n.species,l] %*% P[1:n.species,1:n.species,l] %*% Lambda[1:n.species,1:n.species,l]  
    Q[1:n.species,1:n.species,l] ~ dinvwish(S = I[1:n.species,1:n.species,l], df = n.species + 1)
    P[1:n.species,1:n.species,l] <- Delta[1:n.species,1:n.species,l] %*% Q[1:n.species,1:n.species,l] %*% Delta[1:n.species,1:n.species,l]
    
    for (i in 1:n.species){
      for (j in 1:n.species){
        rho[i,j,l] <- Sigma[i,j,l]/sqrt(Sigma[i,i,l] * Sigma[j,j,l])   
      }
    }
    for(k in 1:n.species){
      for(h in 1:n.year){
        pred1[k,l,h] <- inprod(beta.trend[k,l,1:K], ZZ[h,1:K,k,l])
        mu[k,l,h] <- lbo1[k,l] + inprod(beta.trend[k,l,1:K], ZZ[h,1:K,k,l]) + beta.general[k,l] * res[h] + beta.income[k,l] * REL.COST[h] + beta.drought2[k,l] * wpdsi[h,l] + beta.jobs[k,l] * une[h]
      }
    }
  }
  sig.une ~ dunif(0,5)
  for(h in 1:(n.year)){
    une[h] ~ dnorm(0, sd = sig.une)
  }
  for(l in 1:n.region){
    sig.wpdsi[l] ~ dunif(0,5)
    for(h in 1:n.year){
      wpdsi[h,l] ~ dnorm(0, sd = sig.wpdsi[l])
      hunt.eps[1:n.species,l,h] ~ dmnorm( mu[1:n.species,l,h], cov =  Sigma[1:n.species,1:n.species,l] )
    }
  }
  for (k in 1:n.species){
    for(l in 1:n.region){
      beta.drought2[k,l] ~ dnorm(mu.drought2[l], sd = sig.drought2[l]) 
      beta.jobs[k,l] ~ dnorm(mu.jobs[l], sd = sig.jobs[l]) 
      beta.general[k,l] ~ dnorm(mu.gen[l], sd  = sig.gen[l])
      beta.income[k,l] ~ dnorm(mu.incom[l], sd  = sig.incom[l])
      lbo1[k,l] ~ dnorm(5, sd = 1)
      log(H[k,1,l]) <- hunt.eps[k,l,1]
      z[k,1,l] ~  dpois(H[k,1,l])
      for(h in 2:(n.year)){
        log(H[k,h,l]) <- hunt.eps[k,l,h]
        z[k,h,l] ~  dpois(H[k,h,l])
      }
    }
  }
  for(l in 1:n.region){
    mu.incom[l] ~ dnorm(0, 0.01)
    mu.gen[l] ~ dnorm(0, 0.01)
    mu.drought2[l] ~ dnorm(0, 0.01)
    sig.gen[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.incom[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sig.drought2[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    mu.jobs[l] ~ dnorm(0, 0.01)
    sig.jobs[l] ~ T(dt(0, pow(2.5,-2), 1),0,)
    for (i in 1:n.species){
      sig.trend[i,l] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(k in 1:K){
        
        beta.trend[i,l,k] ~ dnorm(0, sd = sig.trend[i,l])
      }
    } 
    for(h in 2:(n.year)){
      for (i in 1:n.species){
        lambda1[i,h-1,l] <- H[i,h,l]/H[i,h-1,l]
      }
    }
  }
  
})