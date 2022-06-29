#### Modeled Variation in N instead of log.r.harv
#### Spline to deal with serial autocorrelation in N/H
#### No Covariates
#### Added Chukar Site Abundance Back in

### Model Code
code <- nimbleCode( {
  ################################################################################
  ################################################################################
  ### Hunter Effort ###
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
    # alpha.chuk[r] ~ dnorm(0, sd = 100) #Intercept
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