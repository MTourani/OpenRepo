nmodel <- nimbleCode( { 
  
  # Hyperpriors community-level parameters
  mu.lpsi ~ dnorm(0, 0.001)        # Community mean of occupancy (logit)
  mu.lp ~ dnorm(0, 0.001)          # Community mean of detection (logit)
  
  sd.lpsi ~ dunif(0, 10)            # Species heterogeneity in logit(psi)
  sd.lp ~ dunif(0, 10)              # Species heterogeneity in logit(p)                
  
  
  
  # Heterogeneity in detection #################################
  # Community-level fixed effect
  betalp   ~ dnorm(0, 0.1)
  
  # Species-level random effect
  for (k in 1:nspec) {
    
    for (m in 1:ncamar[k]) {
      
      cint.p[m, k] ~ dnorm(0, sd = sd.camar.p)
      
    }
  }
  
  sd.camar.p ~ dunif(0, 10)
  
  
  
  
  # Heterogeneity in occupancy #################################
  # Community-level fixed effects
  betalpsiTrait   ~ dnorm(0, 0.1)
  
  # Interactions
  betalpsiICT   ~ dnorm(0, 0.1)
  betalpsiILT   ~ dnorm(0, 0.1)
  betalpsiICLT   ~ dnorm(0, 0.1)
  
  # Species-level random effects
  for (k in 1:nspec) {
    
    for (m in 1:nhex[k]) {
      
      cint.psi[m, k] ~ dnorm(0, sd = sd.hex.psi)
      
    }
  }
  
  sd.hex.psi ~ dunif(0, 10)
  
  # Species-level fixed effects
  for (k in 1:nspec) {
    
    betalpsi1[k] <- 0
    betalpsiYr[k] ~ dnorm(mu.betalpsiYr, sd = sd.betalpsiYr)
    betalpsiLc[k] ~ dnorm(mu.betalpsiLc, sd = sd.betalpsiLc)
    betalpsiClim[k] ~ dnorm(mu.betalpsiClim, sd = sd.betalpsiClim)
    betalpsiICL[k] ~ dnorm(mu.betalpsiICL, sd = sd.betalpsiICL)
  }
  
  
  # Hyperpriors of heterogeneity in occupancy
  mu.betalpsiYr ~ dnorm(0, 0.001)    # year effect
  sd.betalpsiYr ~ dunif(0, 10)
  
  mu.betalpsiLc ~ dnorm(0, 0.001)   # LC effect
  sd.betalpsiLc ~ dunif(0, 10)
  
  mu.betalpsiClim ~ dnorm(0, 0.001)  # Clim effect
  sd.betalpsiClim ~ dunif(0, 10)
  
  mu.betalpsiICL ~ dnorm(0, 0.001)   # Species-specific interaction effect
  sd.betalpsiICL ~ dunif(0, 10)
  
  
  
  # Priors for species-specific effects in occupancy and detection
  for (k in 1:nspec) {      
    
    lpsi[k] ~ dnorm(mu.lpsi, sd.lpsi)
    lp[k] ~ dnorm(mu.lp, sd.lp)
    
  }
  
  # Ecological submodel for latent occurrence
  for (k in 1:nspec) {
    
    for (i in 1:nsite[k]) {
      
      logit(psi[i, k]) <- lpsi[k] +
        
        # Species-level fixed effects
        betalpsiYr[k] * Xyear[siteID[i, k], 2] +
        betalpsiLc[k] * lc[siteID[i, k], 2] +
        betalpsiClim[k] * clim[siteID[i, k], 3] + 
        betalpsiICL[k] * lc[siteID[i, k], 2] * clim[siteID[i, k], 3] +
        # Community-level fixed effects
        betalpsiTrait * trait[k, 2] + 
        betalpsiICT * clim[siteID[i, k], 3] * trait[k, 2]  +
        betalpsiILT * lc[siteID[i, k], 2] * trait[k, 2] +
        betalpsiICLT * lc[siteID[i, k], 2] * clim[siteID[i, k], 3] * trait[k, 2] +
        # Species-level random effects
        cint.psi[hs[i, k], k]  
      
      
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