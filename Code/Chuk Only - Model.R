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
  
  ################################################################################
  ### Hunter Effort ###
  mu.econ ~ dnorm(0, sd = 10)
  sig.econ ~ T(dt(0, pow(2.5, -2), 1), 0, )
  # beta.econ.hunt ~ dnorm(0, sd = sig.econ)
  
  for(c in 1:n.counties){
    sig.H[c] ~ T(dt(0, pow(2.5,-2), 1),0,)
    alpha.hunt[c] ~ dnorm(0, sd = 5)
    beta.econ.hunt[c] ~ dnorm(mu.econ, sd = sig.econ)
    for(k in 1:K){
      beta.spl.hunt[c,k] ~ dnorm(0, sd = sig.spl.hunt[c])
    } #k
    sig.spl.hunt[c] ~ T(dt(0, pow(2.5,-2), 1),0,)
    
    for(t in 1:n.year){ 
      latent.trend[c,t] <- inprod(beta.spl.hunt[c,1:K], basis[t,1:K]) #spline smoothing
      #Unlinked estimate of Hunter Numbers
      mu.hunt[c,t] <- alpha.hunt[c] + #intercept
        beta.econ.hunt[c] * pred.econ.prime[t] + #SEM economic indicator
        latent.trend[c,t]
      H[c,t] <- exp(hunt.eps[c,t]) #Log Link
      n.hunt[c,t] ~ dnorm(H[c,t], sd = sig.H[c]) #Number of hunters follows Poisson
    } #t
    
    for(t in 1:(n.year-1)){
      lambda.hunt[c,t] <- H[c,t+1]/H[c,t]
      log.r.hunt[c,t] <- log(lambda.hunt[c,t])
    } #t
  }#c
  
  for(t in 1:n.year){
    hunt.eps[1:n.counties,t] ~ dmnorm(mu.hunt[1:n.counties,t], cov =  Sigma.hunt[1:n.counties,1:n.counties] )
  } #t
  
  # Correlation Matrices
  Q.hunt[1:n.counties,1:n.counties] ~ dinvwish(S = I.hunt[1:n.counties,1:n.counties], df = n.counties + 1)
    
  for(c in 1:n.counties){
    sig.hunt[c] ~ dgamma(1,1)
    Delta.hunt[c,c] <- pow(Q.hunt[c,c], -0.5)
    Lambda.hunt[c,c] <- sig.hunt[c]
  } #c
  
  for (c1 in 2:n.counties){
    for (c2 in 1:(c1-1)){
      Lambda.hunt[c1,c2] <- 0
      Delta.hunt[c1,c2] <- 0
    } #c2
  } #c1
  
  Sigma.hunt[1:n.counties,1:n.counties] <- Lambda.hunt[1:n.counties,1:n.counties] %*% P.hunt[1:n.counties,1:n.counties] %*% Lambda.hunt[1:n.counties,1:n.counties]  
  P.hunt[1:n.counties,1:n.counties] <- Delta.hunt[1:n.counties,1:n.counties] %*% Q.hunt[1:n.counties,1:n.counties] %*% Delta.hunt[1:n.counties,1:n.counties]
  
  for (c1 in 1:n.counties){
    for (c2 in 1:n.counties){
      rho.hunt[c1,c2] <- Sigma.hunt[c1,c2]/sqrt(Sigma.hunt[c1,c1] * Sigma.hunt[c2,c2])   
    } #c2
  } #c1
  
  
  ################################################################################
  ### Total Harvest ###
  # mu.wintsev.harv ~ dnorm(0, 0.01)
  # mu.bbs ~ dnorm(0, 0.01)
  mu.hunter.harv ~ dnorm(0, 0.01)
  sig.wintsev.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  sig.bbs ~ T(dt(0, pow(2.5,-2), 1),0,)
  sig.hunter.harv ~ T(dt(0, pow(2.5,-2), 1),0,)
  beta.wintsev.harv ~ dnorm(0, sd  = sig.wintsev.harv)
  beta.bbs.harv ~ dnorm(0, sd  = sig.bbs)
  # beta.hunter.harv ~ dnorm(0, sd = sig.hunter.harv)
  
  for(c in 1:n.counties){
    # beta.wintsev.harv[c] ~ dnorm(mu.wintsev.harv, sd  = sig.wintsev.harv)
    # beta.bbs.harv[c] ~ dnorm(mu.bbs, sd  = sig.bbs)
    beta.hunter.harv[c] ~ dnorm(mu.hunter.harv, sd = sig.hunter.harv)
    sig.N[c] ~ T(dt(0, pow(2.5,-2), 1),0,)
    alpha.harv[c] ~ dnorm(0, sd = 7)
    
    # Process Model
    for(t in 1:n.year){
      #Unlinked estimate of Hunter Numbers
      mu.harv[c,t] <- alpha.harv[c] + # intercepts
        beta.hunter.harv[c] * latent.trend[c,t] + # Latent Hunter Trend
        beta.wintsev.harv * awssi[t,reg.county[c]]  + # Previous winter severity (Affecting Survival)
        beta.bbs.harv * pred.bbs.prime[t]
      
      N[c,t] <- exp(harv.eps[c,t]) #Log Link
      n.harv[c,t] ~ dnorm(N[c,t], sd = sig.N[c]) #Number of hunters follows Poisson
    } #t
    
    for(t in 1:(n.year-1)){
      lambda.harv[c,t] <- N[c,t+1]/N[c,t]
      log.r.harv[c,t] <- log(lambda.harv[c,t])
    } #t-1
  } #c
  
  for(t in 1:n.year){
    harv.eps[1:n.counties,t] ~ dmnorm(mu.harv[1:n.counties,t], cov =  Sigma.harv[1:n.counties,1:n.counties])
  } #t
  
  ### Correlation Matrices
  Q.harv[1:n.counties,1:n.counties] ~ dinvwish(S = I.harv[1:n.counties,1:n.counties], df = n.counties + 1)
  
  for(c in 1:n.counties){
    sig.harv[c] ~ dgamma(1,1)
    Delta.harv[c,c] <- pow(Q.harv[c,c], -0.5)
    Lambda.harv[c,c] <- sig.harv[c]
  } #c
  
  for (c1 in 2:n.counties){
    for (c2 in 1:(c1-1)){
      Lambda.harv[c1,c2] <- 0
      Delta.harv[c1,c2] <- 0
    } #s2
  } #s1
  
  Sigma.harv[1:n.counties,1:n.counties] <- Lambda.harv[1:n.counties,1:n.counties] %*% P.harv[1:n.counties,1:n.counties] %*% Lambda.harv[1:n.counties,1:n.counties]  
  P.harv[1:n.counties,1:n.counties] <- Delta.harv[1:n.counties,1:n.counties] %*% Q.harv[1:n.counties,1:n.counties] %*% Delta.harv[1:n.counties,1:n.counties]  
  
  for (c1 in 1:n.counties){
    for (c2 in 1:n.counties){
      rho.harv[c1,c2] <- Sigma.harv[c1,c2]/sqrt(Sigma.harv[c1,c1] * Sigma.harv[c2,c2])   
    } #c2
  } #c1
  
  
  ################################################################################
  ### Chukar Site Abundance ###
  theta.chuk ~ T(dt(0, pow(2.5,-2), 1),0,) #NB "size" parameter
  
  # for(c in 1:n.counties){
    # mod.chuk[c] ~ dlogis(0,1)
  # }
  
  for(p in 1:n.site){
    C.chuk[p,1] ~ dpois(n.chuk[p,1]) #Equivalent of Poisson lambda
    mod.chuk[p] ~ dlogis(0,1)
    for(t in 1:(n.year.chuk-1)){
      log.r.chuk[p,t] <-  mod.chuk[p] * log.r.harv[county.site[p], t]
      C.chuk[p,t+1] <- exp(log.r.chuk[p,t]) * C.chuk[p,t] #Equivalent of Poisson lambda
    }
     
    for(t in 2:n.year.chuk){ 
      rate.chuk[p,t-1] <- theta.chuk/(theta.chuk + C.chuk[p,t]) #NB success parameter
      n.chuk[p,t] ~ dnegbin(prob = rate.chuk[p,t-1], size = theta.chuk) #obs. # of chukars follow neg-bin
    } #t
  } #p
  
  ################################################################################
  ### Birds per Hunter
  for(t in 1:n.year){
    for(c in 1:n.counties){
      BPH[c,t] <-  N[c,t]/H[c,t]
      BPH2[c,t]<- exp(mu.harv[c,t]-mu.hunt[c,t])
    } #c
  } #t
})