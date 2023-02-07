

nmodel <- nimbleCode( { 
  
  # Hyperpriors community-level parameters
  mu.lpsi ~ dnorm(0, 0.001)        # Community mean of occupancy (logit)
  mu.lp ~ dnorm(0, 0.001)          # Community mean of detection (logit)
  
  sd.lpsi ~ dunif(0, 10)            # Species heterogeneity in logit(psi)
  sd.lp ~ dunif(0, 10)              # Species heterogeneity in logit(p)                
  
  
  # Heterogeneity in detection ######################################
  # Community-level Fixed effects
  betalp   ~ dnorm(0, 0.1)
  
  # Species-level random effects
  for (k in 1:nspec) {
    
    for (m in 1:ncamar[k]) {
      
      cint.p[m, k] ~ dnorm(0, sd = sd.camar.p)
      cint.psi[m, k] ~ dnorm(0, sd = sd.camar.psi)
      
    }
  }
  
  sd.camar.p ~ dunif(0, 10)
  sd.camar.psi ~ dunif(0, 10)
  
  

  # Heterogeneity in occupancy ######################################
  
  # Species level random effects
  for (k in 1:nspec) {
    
    for (m in 1:nhex[k]) {
      
      cint.betaLc[m, k] ~ dnorm(0, sd =  sd.betalpsiLc[k])
      
      
    }
  }
  

  mu.mubeta ~ dnorm(0, 0.001)
  sd.mubeta ~ dunif(0, 10)
  a.sdbeta  ~ dgamma(0.5, 0.0005)
  b.sdbeta ~ dgamma(0.5, 0.0005)
  
  
  for (k in 1:nspec) {
    
    betalpsi1[k] <- 0
    betalpsiYr[k] ~ dnorm(mu.betalpsiYr, sd = sd.betalpsiYr)

    mu.betalpsiLc[k] ~ dnorm(mu.mubeta, sd = sd.mubeta)  
    sd.betalpsiLc[k] ~ dgamma(a.sdbeta, b.sdbeta)
    
    
  }
  
  
  # Hyperpriors of heterogeneity in occupancy
  mu.betalpsiYr ~ dnorm(0, 0.001) # year effect
  sd.betalpsiYr ~ dunif(0, 10)
  
  
  # Priors for species-specific effects in occupancy and detection
  for (k in 1:nspec) {      
    
    lpsi[k] ~ dnorm(mu.lpsi, sd.lpsi)
    lp[k] ~ dnorm(mu.lp, sd.lp)
    
  }
  
  # Ecological submodel for latent occurrence
  for (k in 1:nspec) {
    
    for (i in 1:nsite[k]) {
      
      # Species level hex-specific covariate effects
      betaLc[i, k] <- mu.betalpsiLc[k] + cint.betaLc[hs[i, k], k]
      
      logit(psi[i, k]) <- lpsi[k] + betalpsiYr[k] * Xyear[siteID[i, k], 2] 
                                  + betaLc[i, k] * Lc[siteID[i, k], 1] 
                                  + cint.psi[ms[i, k], k] 
      
      
      z[i, k] ~ dbern(psi[i, k])
      
    }
    
  }
  
  # Observation submodel
  for (k in 1:nspec) { 
    
    for (i in 1:nsite[k]) {
      
      logit(p[i, k]) <- lp[k] +  betalp * acam[siteID[i, k]] + cint.p[ms[i, k], k] 
      
      y[i, k] ~ dbin(z[i, k] * p[i,k], nrep[i])
      
    }
    
  }
  
  
  
})