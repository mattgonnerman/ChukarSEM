### Model Code
code <- nimbleCode( {
  ################################################################################
  ### Predictors
  #Personal Disposable Income
  alpha.pdi ~ dnorm(0, 0.001) #intercept
  beta.t.pdi ~ dnorm(0, 0.01) #year
  sig.pdi~ T(dt(0, pow(2.5,-2), 1),0,)
  
  ar1.pdi ~ dunif(-1,1) #Autoregressive parameter
  
  pdi.trend[1] <- alpha.pdi + beta.t.pdi * 1
  mu.pdi[1] <- pdi.trend[1]
  for(t in 2:n.year){
    pdi.trend[t] <- alpha.pdi + beta.t.pdi * t
    mu.pdi[t] <- pdi.trend[t] + ar1.pdi * (PDI[t-1] - pdi.trend[t-1])
  } #t
  
  #Gas Prices
  alpha.gas ~ dnorm(1.5, 1)
  sig.gas~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.gas[1] ~ dnorm(0, 0.01)
  beta.gas[2] ~ dnorm(0, 0.01)
  
  #Relative Cost of Gas
  for(t in 1:n.year){
    PDI[t] ~ dnorm(mu.pdi[t], sd = sig.pdi)
    GAS[t] ~ T(dnorm(alpha.gas + beta.gas[era.gas[t]]*t, sd = sig.gas),0,)
    rel.cost[t] <- ((PDI[t]/GAS[t]) - 2.581635)/0.8894599
  } #t
  
  #Unemployment Rate
  sig.une ~ dunif(0,5)
  for(t in 1:(n.year)){
    une[t] ~ dnorm(0, sd = sig.une)
  } #t
  
  #Resident License sales
  sig.res ~ dunif(0,5)
  for(t in 1:n.year){
    res[t] ~ dnorm(0, sd = sig.res)
  }
  
  #Drought Index
  for(r in 1:n.region){
    sig.wpdsi[r] ~ dunif(0,5)
    sig.pdsi[r] ~ dunif(0,5)
    for(t in 1:n.year){
      wpdsi[t,r] ~ dnorm(0, sd = sig.wpdsi[r])
      pdsi[t,r] ~ dnorm(0, sd = sig.pdsi[r])
    } #t
  } #r
  
  #Winter Severity
  sig.awssi~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.awssi[1] ~ dnorm(0, 0.01)
  beta.awssi[2] ~ dnorm(0, 0.01)
  for(r in 1:n.region){
    alpha.awssi[r] ~ dnorm(1.5, 1)
    for(t in 1:n.year){
      awssi[r,t] ~ dnorm(alpha.awssi[r] + beta.awssi[era.awssi[t]]*t, sd = sig.awssi)
    }
  }
  
  #Rabbits
  sig.rabbits~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.rabbits ~ dnorm(0, 0.01)
  for(r in 1:n.region){
    alpha.rabbits[r] ~ dnorm(1.5, 1)
    for(t in 1:n.year){
      rabbits[t,r] ~ dnorm(alpha.rabbits[r] + beta.rabbits*t, sd = sig.rabbits)
    }
  }
  
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
  
  #Northern Harrier
  sig.nhar ~ dunif(0,5)
  for(t in 1:(n.year)){
    nharrier[t] ~ dnorm(0, sd = sig.nhar)
  } #t
  
  ################################################################################
  ### Hunter Effort ###
  mu.drought.hunt ~ dnorm(0, 0.01)
  sig.drought.hunt ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.wintsev.hunt ~ dnorm(0, 0.01)
  sig.wintsev.hunt ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.jobs ~ dnorm(0, 0.01)
  sig.jobs ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.income ~ dnorm(0, 0.01)
  sig.income ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.license ~ dnorm(0, 0.01)
  sig.license ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
    beta.drought.hunt[s] ~ dnorm(mu.drought.hunt, sd = sig.drought.hunt)
    beta.wintsev.hunt[s] ~ dnorm(mu.wintsev.hunt, sd  = sig.wintsev.hunt)
    beta.jobs[s] ~ dnorm(mu.jobs, sd = sig.jobs)
    beta.income[s] ~ dnorm(mu.income, sd  = sig.income)
    beta.license[s] ~ dnorm(mu.license, sd  = sig.license)
  }
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.hunt[s,r] ~ dnorm(5, sd = 1)
      for(k in 1:K){
        beta.spl.hunt[s,r,k] ~ dnorm(0, sd = sig.spl.hunt[s,r])
      } #k
      sig.spl.hunt[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      for(t in 1:n.year){ 
        #Unlinked estimate of Hunter Numbers
        mu.hunt[s,r,t] <- alpha.hunt[s,r] + #intercept
          beta.drought.hunt[s] * wpdsi[t,r] + #concurrent winter drought index
          beta.wintsev.hunt[s] * awssi[r,t] + #concurrent winter severity
          beta.jobs[s] * une[t] + #concurrent years unemployment
          beta.income[s] * rel.cost[t] + #PDI/Gas Price
          beta.license[s] * res[t] + #Hunting licences sold that season
          inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #spline smoothing
        
        pred.spl.hunt[s,r,t] <- inprod(beta.spl.hunt[s,r,1:K], Z.hunt[t,1:K,s,r]) #Derive spline smoothing for examination later
        
        H[s,t,r] <- exp(hunt.eps[s,r,t]) #Log Link
        
        n.hunt[s,t,r] ~ dpois(H[s,t,r]) #Number of hunters follows Poisson
      } #t
      
      for(t in 1:(n.year-1)){
        lambda.hunt[s,t,r] <- H[s,t+1,r]/H[s,t,r]
        log.r.hunt[s,t,r] <- log(lambda.hunt[s,t,r])
      }
    } #s
    
    for(t in 1:n.year){
      hunt.eps[1:n.species,r,t] ~ dmnorm(mu.hunt[1:n.species,r,t], cov =  Sigma.hunt[1:n.species,1:n.species,r] )
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
  mu.drought.harv ~ dnorm(0, 0.01)
  sig.drought.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.wintsev.harv ~ dnorm(0, 0.01)
  sig.wintsev.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.rabbit ~ dnorm(0, 0.01)
  sig.rabbit ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.raven ~ dnorm(0, 0.01)
  sig.raven ~ T(dt(0, pow(2.5,-2), 1),0,)
  mu.nharrier ~ dnorm(0, 0.01)
  sig.nharrier ~ T(dt(0, pow(2.5,-2), 1),0,)
  
  for(s in 1:n.species){
    beta.drought.harv[s] ~ dnorm(mu.drought.harv, sd = sig.drought.harv)
    beta.wintsev.harv[s] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
    beta.rabbit.harv[s] ~ dnorm(mu.rabbit, sd  = sig.rabbit)
    beta.raven.harv[s] ~ dnorm(mu.raven, sd  = sig.raven)
    beta.nharrier.harv[s] ~ dnorm(mu.nharrier, sd  = sig.nharrier)
  } #s
  
  for(r in 1:n.region){
    for(s in 1:n.species){
      alpha.harv[s,r] ~ dnorm(0, sd = 1)
      for(k in 1:K){
        beta.spl.harv[s,r,k] ~ dnorm(0, sd = sig.spl.harv[s,r])
      } #k
      sig.spl.harv[s,r] ~ T(dt(0, pow(2.5,-2), 1),0,)
      
      # Process Model
      N[s,1,r] ~ dpois(n.harv[s,1,r]) #Total harvest, Year 1
      
      for(t in 2:(n.year)){
        mu.harv[s,t-1,r] <- alpha.harv[s,r] +#regression formula
          beta.drought.harv[s] * pdsi[t-1,r] + #previous breeding season drought index
          beta.wintsev.harv[s] * awssi[r,t-1] + #concurrent winter severity
          beta.rabbit.harv[s] * rabbits[t-1,r] + #concurrent rabbit harvest
          beta.raven.harv[s] * raven[t] + #prior BBS index 
          beta.nharrier.harv[s] * nharrier[t] + #prior BBS index 
          inprod(beta.spl.harv[s,r,1:K], Z.harv[t,1:K,s,r]) #spline smoothing
        
        pred.spl.harv[s,r,t-1] <- inprod(beta.spl.harv[s,r,1:K], Z.harv[t,1:K,s,r]) #Derive spline smoothing for examination later
        
        lambda.harv[s,t-1,r] <- exp(log.r.harv[s,t-1,r]) #link function
        N[s,t,r] <- lambda.harv[s,t-1,r] * N[s,t-1,r] #number available = change since last year
        n.harv[s,t,r] ~  dpois(N[s,t,r]) #Number harvested follows Poisson
      } #t
    } #s 
    
    for(t in 2:(n.year)){
      #Change in total harvest, log.r.harv[t=1] is 1976-1977
      log.r.harv[1:n.species,t-1,r]  ~ dmnorm(mu.harv[1:n.species,t-1,r],
                                              cov =  Sigma.harv[1:n.species,1:n.species,r])
    } #t
    
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
  ### Sage Grouse Wing-Bee ###
  beta.drought.sg ~ dnorm(0, 0.01)
  beta.wintsev.sg ~ dnorm(0, 0.01)
  beta.rabbit.sg ~ dnorm(0, 0.01)
  beta.raven.sg ~ dnorm(0, 0.01)
  beta.nharrier.sg ~ dnorm(0, 0.01)
  
  for(r in 1:n.region){
    alpha.sg[r] ~ dnorm(0, sd = 100) #Intercept
    mod.sg[r] ~ dlogis(0,1) #Constant modifier to translate harv change to wingb change
    theta.sg[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    
    for(t in 1:n.years.sg){
      sg.eps[r, t] <- mod.sg[r] * log.r.harv[7, t+27, r] #log.r.harv[t=25] is 2004
      
      log.r.sg[r,t] <- alpha.sg[r] + 
        sg.eps[r,t] + #Unlinked change in recruitment
        beta.drought.sg * pdsi[t+28,r] + #previous breeding season drought index
        beta.wintsev.sg * awssi[r,t+27] + #previous winter severity
        beta.rabbit.sg * rabbits[t+27,r] + #previous year's rabbit harvest
        beta.raven.sg * raven[t+29] + #previous spring BBS raven index
        beta.nharrier.sg * nharrier[t+29]  #previous spring BBS northern harrier index
      
      rate.sg[r,t] <- theta.sg[r]/(theta.sg[r] + (AHY.sg[r,t]*exp(log.r.sg[r,t]))) #NB rate
      HY.sg[r,t] ~ dnegbin(prob = rate.sg[r,t], size = theta.sg[r])
    } #t
  } #r
  
  ################################################################################
  ### Chukar Site Abundance ###
  beta.drought.chuk ~ dnorm(0, 0.01)
  beta.wintsev.chuk ~ dnorm(0, 0.01)
  beta.rabbit.chuk ~ dnorm(0, 0.01)
  beta.raven.chuk ~ dnorm(0, 0.01)
  beta.nharrier.chuk ~ dnorm(0, 0.01)
  
  for(r in 1:n.region){
    alpha.chuk[r] ~ dnorm(0, sd = 100) #Intercept
    theta.chuk[r] ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
    mod.chuk[r] ~ dlogis(0,1)
  }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    
    for(t in 2:n.year.chuk){
      chuk.eps[p, t-1] <- mod.chuk[reg.chuk[p]] * log.r.harv[3, t+13, reg.chuk[p]] #log.r.harv[t=13] is 1990-1991
      
      log.r.chuk[p,t-1] <- alpha.chuk[reg.chuk[p]] + 
        chuk.eps[p,t-1] + #Unlinked change in abundance, log.harv.chuk[t=1] is 1990-1991
        beta.drought.chuk * pdsi[t+12,reg.chuk[p]] + #previous breeding season drought index
        beta.wintsev.chuk * awssi[reg.chuk[p],t+12] + #previous year's winter severity
        beta.rabbit.chuk * rabbits[t+12,reg.chuk[p]] + #previous year's rabbit harvest
        beta.raven.chuk * raven[t+13] + #previous year's spring BBS index
        beta.nharrier.chuk * nharrier[t+13] #previous year's spring BBS index
      
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