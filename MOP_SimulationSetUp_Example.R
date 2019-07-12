rm(list=ls())
gc()
## ------ IMPORT REQUIRED LIBRARIES ------
.libPaths(c("C:/Rfolder", .libPaths())) 
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(nimble)
library(rjags)
library(spdep)
library(maptools)
library(stringr)
library(abind)
library(R.utils)
## ------ SET REQUIRED WORKING DIRECTORIES ------
source("C:/Users/admin/Documents/PROJECTS/RovQuant/Temp/MT/myWorkingDirectories.R")

## ------ SOURCE THE REQUIRED FUNCTIONS ------

sourceDirectory("C:/Users/admin/Documents/PROJECTS/RovQuant/Source",modifiedOnly=FALSE) 
sourceDirectory("C:/Users/admin/Documents/PROJECTS/RovQuant/Source_Nimble",modifiedOnly=FALSE) 
source("C:/Users/admin/Documents/PROJECTS/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/makeSxyInits.R")
source("C:/Users/admin/Documents/PROJECTS/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/calculateDetProb_a0.R")
source("C:/Users/admin/Documents/PROJECTS/RovQuant/Temp/MT/HR/FUNCTIONS/dbin_LESSCachedSparsea0.R")

## ---------------------------------------------------------------------------------------------
## ------ I.LOAD RAW DATA ------
dir.dropbox <- "C:/Users/admin/Dropbox (AQEG)/AQEG Team Folder/RovQuant/"
setwd(dir.dropbox) 

## ------ SOURCE THE REQUIRED FUNCTIONS ------
sourceDirectory(dir.function, modifiedOnly = FALSE)
sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)

inv.logit <- function(inValues){1.0 / (1.0 + exp(-inValues))}

## -----------------------------------------------------------------------------------------------------------
## ------ I. SIMULATE DATA ------

mySim <- list( HABITAT = list( resolution = 1#parm.df2$habitat.resolution[i]
                               , buffer = 0#parm.df2$buffer[i]
                               , grid.size = 50#parm.df2$grid.size[i])
)
               , DETECTORS = list( resolutions = 1.2#parm.df2$detector.resolutions[i]
                                   , divisions = 1#parm.df2$detector.divisions[i]
                                   , resolutions.cam = 1#parm.df2$detector.resolutions.cam[i])
               )
               , DETECTIONS = list( N = 30#parm.df2$N[i]
                                    , a0 = 1#parm.df2$a0[i]
                                    , a0.cam = 0.75#parm.df2$a0.cam[i]
                                    , sigma = 1#parm.df2$sigma[i]
                                    , aug.factor = 1#parm.df2$aug.factor[i]
                                    , a = 0.25#parm.df2$a[i])
                                    )
)


### ==== 1.GENERATE HABITAT ====
myStudyArea <- MakePolygon(max.x = mySim$HABITAT$grid.size, max.y = mySim$HABITAT$grid.size, random = FALSE)
myHabitat <- MakeHabitat(poly = myStudyArea, resolution = mySim$HABITAT$resolution, buffer = mySim$HABITAT$buffer)                                 

### ==== 2.GENERATE DETECTORS ====
myDetectors <- MakeSearchGrid(data = myStudyArea, resolution = mySim$DETECTORS$resolutions)  
myDetectors.cam <- MakeSearchGrid(data = myStudyArea, resolution = mySim$DETECTORS$resolutions.cam)  

### ==== 3.RESCALE DETECTORS & HABITAT ====
scaled <- UTMToGrid(data.sp = myDetectors$detector.sp, grid.sp = myHabitat$habitat.sp, plot.check = F)
scaled.cam <- UTMToGrid(data.sp = myDetectors.cam$detector.sp, grid.sp = myHabitat$habitat.sp, plot.check = F)

### ==== 4.SIMULATE AC LOCATIONS ====
mySimulatedACs <- SimulateACs(N = mySim$DETECTIONS$N, habitat.sp = myHabitat$habitat.sp, plot = TRUE)

### ==== 5.SIMULATE DETECTIONS ====
myDetections <- SimulateDetection( a0 = mySim$DETECTIONS$a0
                                   , sigma = mySim$DETECTIONS$sigma
                                   , AC.sp = mySimulatedACs
                                   , detector.sp = myDetectors$detector.sp
                                   , type = "Bernoulli"
                                   , plot = TRUE)

myDetections.cam <- SimulateDetection( a0 = mySim$DETECTIONS$a0.cam
                                       , sigma = mySim$DETECTIONS$sigma
                                       , AC.sp = mySimulatedACs
                                       , detector.sp = myDetectors.cam$detector.sp
                                       , type = "Bernoulli"
                                       , plot = TRUE)

### ==== 6.GENERATE THE DIFFERENT TYPES OF DETECTIONS ====
y.ORIGINAL.ngs <- myDetections$y.all

identifiedDets <- apply(y.ORIGINAL.ngs, c(1,2), function(x)rbinom(1, 1, mySim$DETECTIONS$a))

y.IDENTIFIED.ngs <- ifelse(y.ORIGINAL.ngs == 1 & identifiedDets == 1, 1, 0)

y.POOLED.ngs <- as.numeric(apply(y.ORIGINAL.ngs, 2, function(x)any(x >= 1)))

y.UNIDENTIFIED.ngs <- ifelse(y.ORIGINAL.ngs == 1 & identifiedDets == 0, 1, 0)
y.UNIDENTIFIED.ngs <- as.numeric(apply(y.UNIDENTIFIED.ngs, 2, function(x)any(x >= 1)))  

y.cam <- apply(myDetections.cam$y.all, 2, function(x)any(x >= 1))

### ==== 7.AUGMENT DATA ====
y <- MakeAugmentation(y = y.IDENTIFIED.ngs, aug.factor =1)

### ==== 8.RECONSTRUCT z ====
z <- z.init <- rep(NA, dim(y)[1])
detected <- apply(y, 1, function(x)any(x >= 1))
z[detected] <- 1
z.init[!detected] <- 1

## -------------------------------------------------------------------------------------------------------------
## ----- II. MODEL DEFINITION -----  
## -------------------------------------------------------------------------------------------------------------
MOP <- nimbleCode({
  ##-----------------------------------------------------------------------------------------------
  ##---------- SPATIAL PROCESS ----------## 
  for(i in 1:M){
    sxy[i,1] ~ dunif(0, x.max)
    sxy[i,2] ~ dunif(0, y.max)
  }#i
  
  ##-----------------------------------------------------------------------------------------------
  ##---------- DEMOGRAPHIC PROCESS ----------## 
  psi ~ dunif(0,1)
  for(i in 1:M){
    z[i] ~ dbern(psi)
  }#i
  
  ##-----------------------------------------------------------------------------------------------
  ##---------- DETECTION PROCESS ----------## 
  sigma ~ dunif(0,100)
  p0 ~ dunif(0,1)
  p0.cam ~ dunif(0,1)
  a ~ dunif(0,1)  
  
  ## MARKED NGS
  for(i in 1:M){
    p[i,1:J] <- calculateDetProb(sxy[i,1:2], detector.xy[1:J,1:2], p0, sigma)
    p.id[i,1:J] <- p[i,1:J]*a
    p.unid[i,1:J] <- p[i,1:J]*(1-a)*z[i]
    y[i,1:J] ~ dbern_vector(p.id[i,1:J], z[i])
  }#i
  
  ## UNMARKED NGS 
  for(j in 1:J){
    bigP[j] <- 1 - prod((1-p.unid[1:M,j]))
  }#j
  unid.y.ngs[1:J] ~ dbern_vector(bigP[1:J], 1)
  
  ## CAMERA-TRAPS
  for(i in 1:M){
    p.cam[i,1:J.cam] <- calculateDetProb(sxy[i,1:2], detector.xy.cam[1:J.cam,1:2], p0.cam, sigma)
    p.prim.cam[i,1:J.cam] <- p.cam[i,1:J.cam] * z[i]
  }#i
  
  for(j in 1:J.cam){
    bigP.cam[j] <- 1 - prod((1-p.prim.cam[1:M,j]))
  }#j
  pooled.y.cam[1:J.cam] ~ dbern_vector(bigP.cam[1:J.cam], 1)
  
  ##---------- DERIVED PARAMETERS ----------##
  N <- sum(z[1:M])
})

## -------------------------------------------------------------------------------------------------------------
## ----- III. CREATE INPUT -----

nimData <- list( z = z                                                     
              , y = y                                                     
              , unid.y.ngs = y.UNIDENTIFIED.ngs 
              , pooled.y.cam = y.cam 
              , detector.xy = scaled$data.scaled.xy                       
              , detector.xy.cam = scaled.cam$data.scaled.xy                   
              , x.max = ncol(myHabitat$habitat.mx)
              , y.max = nrow(myHabitat$habitat.mx))

nimConstants <- list( M = dim(y)[1]
                 , J = dim(y)[2]
                 , J.cam = dim(scaled.cam$data.scaled.xy)[1])

## -------------------------------------------------------------------------------------------------------------
### ==== 1.PARAMETERS TO SAVE  ====

params <- c( "N"
           , "sigma"
           , "p0"
           , "p0.cam"
           , "a"
           , "psi"
)


### ==== 2.INITIAL VALUES  ====

inits <- list( "z" = z.init
             , "sigma" = runif(1,0,10) 
             , "p0" = runif(1,0,1)
             , "p0.cam" = runif(1,0,1)
             , "psi" = runif(1,0,1)
             , "a" = runif(1,0,1)
)

## -------------------------------------------------------------------------------------------------------------
## ----- IV. MODEL RUN -----

### ==== 1.MODEL BUILDING  ====
nimble.model <- nimbleModel(code = MOP
                          , constants = nimConstants
                          , data = nimData
                          , inits = inits
                          , check = FALSE)
  
### ==== 2.MODEL COMPILATION  ====
Sys.setenv(PATH = paste("C:/Rfolder/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rfolder/Rtools/mingw_$(WIN)/bin/")
Cmodel <- compileNimble(nimble.model)

### ==== 3.MCMC CONFIGURATION ====
mcmcConf <- configureMCMC(nimble.model, monitors = params)
MCMC <- buildMCMC(mcmcConf)

### ==== 4.MCMC COMPLIATION ====
cMCMC <- compileNimble(MCMC, project = nimble.model)

### ==== 5.MCMC SAMPLING ====
samplesList <- runMCMC(cMCMC
                      , niter = 1000
                      , nburnin = 100
                      , nchains = 3
                      , samplesAsCodaMCMC = TRUE)

## -------------------------------------------------------------------------------------------------------------
## ----- V.OUTPUT -----
summary(mcmc.list(samplesList))


  