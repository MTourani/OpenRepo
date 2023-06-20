
library(nimble)
library(coda)

nmodel <- nimbleCode( { 
  
  # Hyperpriors community-level parameters
  mu.lpsi ~ dnorm(0, 0.001)        # Community mean of occupancy (logit)
  mu.lp ~ dnorm(0, 0.001)          # Community mean of detection (logit)
  
  sd.lpsi ~ dunif(0, 10)            # Species heterogeneity in logit(psi)
  sd.lp ~ dunif(0, 10)              # Species heterogeneity in logit(p)                
  
  
  # Heterogeneity in detection #########################################
  # Community-level Fixed effects
  betalp   ~ dnorm(0, 0.1)
  
  # Species-level random effects
  for (k in 1:nspec) {
    
    for (m in 1:ncamar[k]) {
      
      cint.p[m, k] ~ dnorm(0, sd = sd.camar.p)
      
    }
  }
  
  sd.camar.p ~ dunif(0, 10)
  
  
  
  # Heterogeneity in occupancy #########################################
  
  # Species-level random effects
  for (k in 1:nspec) {
    
    for (m in 1:nhex[k]) {
      
      cint.psi[m, k] ~ dnorm(0, sd = sd.hex.psi)
      
      
    }
  }
  
  sd.hex.psi ~ dunif(0, 10)
  
  
  # Heterogeneity in occupancy
  for (k in 1:nspec) {
    
    betalpsi1[k] <- 0
    betalpsiYr[k] ~ dnorm(mu.betalpsiYr, sd = sd.betalpsiYr)
    betalpsiLc[k] ~ dnorm(mu.betalpsiLc, sd = sd.betalpsiLc)
    betalpsiClim[k] ~ dnorm(mu.betalpsiClim, sd = sd.betalpsiClim)
    betalpsiI[k] ~ dnorm(mu.betalpsiI, sd = sd.betalpsiI)
    
  }
  
  
  # Hyperpriors of heterogeneity in occupancy
  mu.betalpsiYr ~ dnorm(0, 0.001) # year effect
  sd.betalpsiYr ~ dunif(0, 10)
  
  mu.betalpsiLc ~ dnorm(0, 0.001) # For LC effect
  sd.betalpsiLc ~ dunif(0, 10)
  
  mu.betalpsiClim ~ dnorm(0, 0.001) # Clim effect
  sd.betalpsiClim ~ dunif(0, 10)
  
  mu.betalpsiI ~ dnorm(0, 0.001) # Interaction effect
  sd.betalpsiI ~ dunif(0, 10)
  
  
  # Priors for species-specific effects in occupancy and detection
  for (k in 1:nspec) {      
    
    lpsi[k] ~ dnorm(mu.lpsi, sd.lpsi)
    lp[k] ~ dnorm(mu.lp, sd.lp)
    
  }
  
  # Ecological submodel for latent occurrence
  for (k in 1:nspec) {
    
    for (i in 1:nsite[k]) {
      
      
      logit(psi[i, k]) <- lpsi[k] + betalpsiYr[k] * Xyear[siteID[i, k], 2] 
                                  + betalpsiLc[k] * Lc[siteID[i, k], 1] 
                                  + betalpsiClim[k] * clim[siteID[i, k], 3] 
                                  + betalpsiI[k]*Lc[siteID[i, k], 1]*clim[siteID[i, k], 3] 
                                  + cint.psi[hs[i, k], k]  
      
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

### ==== 1.MODEL BUILDING  ====
model <- nimbleModel(      code      = nmodel
                         , constants = constants
                         , data      = ndata
                         , inits     = ninits
                         , check     = FALSE
                         , calculate = FALSE)


### ==== 2.MODEL COMPILATION  ====
Cmodel <- compileNimble(model)

### ==== 3.MCMC CONFIGURATION ====
mcmcConf <- configureMCMC(model, monitors = params)
MCMC <- buildMCMC(mcmcConf)

### ==== 4.MCMC COMPLIATION ====
cMCMC <- compileNimble(MCMC, project = model)

### ==== 5.MCMC SAMPLING ====
samplesList <- runMCMC(    cMCMC
                         , niter   = 1000
                         , nburnin = 1
                         , nchains = 1
                         , samplesAsCodaMCMC = TRUE)
psamples <- do.call(rbind, samplesList)

### ==== 6.Simulate new values from the posterior ====
sim <- simulateFromPosterior( model = model, nodes = dimnames(psamples)[[2]])
csim <- compileNimble(sim, project = cmodel$compiledModel, resetFunctions = TRUE)
first <- dim(samplesList[[1]])[1] - (dim(samplesList[[1]])[1] * 0.5) + 1
last  <- dim(samplesList[[1]])[1]
# Because using all samples will be slow and requires lots of memory
pnout <- do.call(rbind, lapply(samplesList, function(x) x[first:last,]))
csim$run(pnout)

### ==== 7.Calculate discrepancy between observed and expected values ====
bp <- 1:dim(ndata$y)[2]
e <- 0.0001
fit_y         <- matrix(NA, length(csim$mv[["y"]][[1]]), dim(ndata$y)[2])
fit_y_new     <- matrix(NA, length(csim$mv[["y"]][[1]]), dim(ndata$y)[2])
fit_agg_y     <- array(NA, dim = c(length(csim$mv[["y"]][[1]]), dim(ndata$y)[2], dim(ndata$y)[1]))
fit_agg_y_new <- array(NA, dim = c(length(csim$mv[["y"]][[1]]), dim(ndata$y)[2], dim(ndata$y)[1]))

for (i in 1:dim(ndata$y)[2]) {
  for (j in 1:length(csim$mv[["y"]][[1]])) {
    E_agg <- sum(csim$mv[["p"]][[j]][, i, drop = FALSE] * csim$mv[["z"]][[j]][1:max(constants$nsite), i, drop = FALSE], na.rm = TRUE)
    fit_agg_y[j, i, ]     <- (ndata$y[,i] - E_agg)^2 / (E_agg + e)
    fit_agg_y_new[j, i, ] <- (csim$mv[["y"]][[j]][, i] - E_agg)^2 / (E_agg + e)
    fit_y[j, i]       <- sum(fit_agg_y[j, i, ])
    fit_y_new[j, i]   <- sum(fit_agg_y_new[j, i, ])
  }
}

for (k in 1:dim(ndata$y)[2]) {
  bp[k] <- mean(fit_y[,k] > fit_y_new[,k])
}

BP <- mean(apply(fit_y, 1, sum) > apply(fit_y_new, 1, sum))


