


nimModel <- nimbleCode({ 
  ##---------- SPATIAL PROCESS ----------##
  betaBNF ~ dnorm(0.0,0.01) 
  betaNEU ~ dnorm(0.0,0.01) 
  betaBNFp ~ dnorm(0.0,0.01) 
  betaSUMp ~ dnorm(0.0,0.01) 
  betaE ~ dnorm(0.0,0.01)
  betaER ~ dnorm(0.0,0.01)
  betaR ~ dnorm(0.0,0.01)
  
  for(h in 1:n.cells) {                                        #Equation (1)
    habIntensity[h] <- exp(   betaBNF * REGION_BNF[h] +
                                betaBNFp * REGION_BNFp[h] +                             
                                betaNEU * REGION_NEU[h] +
                                betaSUMp * REGION_SUMp[h] +
                                betaR * R[h] +
                                betaE * E[h] +
                                betaER * E[h] * R[h])
  }#h
  
  sumHabIntensity <- sum(habIntensity[1:n.cells])
  logHabIntensity[1:n.cells] <- log(habIntensity[1:n.cells])
  logSumHabIntensity <- log(sumHabIntensity)  
  
  for(i in 1:M) {
    sxy[i,1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.cells, 1:2],
      upperCoords = upperHabCoords[1:n.cells, 1:2],
      logIntensities = logHabIntensity[1:n.cells],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max, 1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  ##---------- DEMOGRAPHIC PROCESS ----------##
  N <- sum(z[1:M])                                          #Equation (2) and (3)
  
  psi ~ dunif(0, 1)                                            
  
  for (i in 1:M) {
    z[i] ~ dbern(psi)
  }#i
  
  ##---------- OBSERVATION PROCESS ----------##
  for(u in 1:n.admin) {                                     #Equation (4) and (5)
    p0[u] ~ dunif(0, 1)
  }#u
  
  betaL ~ dnorm(0, 0.0001)
  
  sigma ~ dunif(0, 30)
  
  for (i in 1:M) {
    y[i, 1:nMaxDetectors] ~  dbinom_sparseLocalSCRTrapCovInt( detNums = nbDetections[i],
                                                              detIndices = yDets[i, 1:nMaxDetectors],
                                                              size = trials[1:J],
                                                              p0 = p0[1:n.admin],
                                                              sigma = sigma,
                                                              s = sxy[i, 1:2],
                                                              trapCoords = detector.xy[1:J, 1:2],
                                                              localTrapsIndices = detectorIndex[1:n.cellsSparse, 1:maxNBDets],
                                                              localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
                                                              resizeFactor = resizeFactor,
                                                              habitatGrid = habitatIDDet[1:y.maxDet, 1:x.maxDet],
                                                              trapCov = DetsCovs[1:n.detectors],
                                                              trapIntIndex = trapIntIndex[1:n.detectors],
                                                              betaTrap = betaL,
                                                              indicator = z[i])
  }#i
  
})
