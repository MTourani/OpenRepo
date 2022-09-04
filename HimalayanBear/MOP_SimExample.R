
## ------ IMPORT REQUIRED LIBRARIES ------

library(nimble)                    
library(R.utils)                   

## ------ Source FUNCTIONS ------
sourceDirectory("~/OpenRepo/FUN/", modifiedOnly = FALSE) 

## ------ DEFINE INPUT DIRECTORY ------
path.in <- "~/OpenRepo/NIMIN/"
path.proc <- "~/OpenRepo/NIMPROC/"
path.out <- "~/OpenRepo/NIMOUT/"

file.list <- list.files(path.in)

## ------------------------------------------------------------------------------------------------------------------ I.LOAD DATA -----  

  ##--Loop over input files
while(length(file.list)>0){
  
  ##--Sample one file at random
  set <- sample(file.list, 1)
  
  ##--Load the sampled data
  load(paste(path.in, set, sep=""))
  print(set)
  file.rename(from = paste(path.in, set, sep = ""), to = paste(path.proc, set, sep = "")) 
  
## ------------------------------------------------------------------------------------------------------------------ II.MODEL RUNS ----- 
  
  ### ==== 1.MODEL BUILDING  ====
  nimble.model <- nimbleModel(  code = mod
                              , constants = constants
                              , data = data
                              , inits = inits
                              , check = FALSE)
  
  
  ### ==== 2.MODEL COMPILATION  ====
  Cmodel <- compileNimble(nimble.model)
  
  ### ==== 3.MCMC CONFIGURATION ====
  mcmcConf <- configureMCMC(nimble.model, monitors = params.mod)
  MCMC <- buildMCMC(mcmcConf)
  
  ### ==== 4.MCMC COMPLIATION ====
  cMCMC <- compileNimble(MCMC, project = nimble.model)
  
  ### ==== 5.MCMC SAMPLING ====
  samplesList <- runMCMC(cMCMC
                         , niter = 10000
                         , nburnin = 1
                         , nchains = 3
                         , samplesAsCodaMCMC = TRUE)
  
  outname <- paste(path.out,"NimOutFor",set, sep = "")
  save(samplesList, file = outname)
  
  file.list <- list.files(path.in)
}

## ------------------------------------------------------------------------------------------------------------------ III.PROCESS THE OUTPUT -----  
par(mfrow=c(1,2))

col.list <- adjustcolor(rev(c("darkblue","lightblue","turquoise","red")), alpha = 0.8)

for(this.mod in 1:4){
  load(paste("~/PROJECTS/OpenRepo/NIMOUT/NimOutForSimDataSet1Rep4mod", this.mod, ".RData", sep = ""))
  if(this.mod == 1)plot(density(cbind(samplesList$chain1[,"N"], samplesList$chain2[,"N"], samplesList$chain3[,"N"])), main ="Posterior N", xlim = c(1,60), ylim = c(0,0.5))
  if(this.mod!=1)lines(density(cbind(samplesList$chain1[,"N"], samplesList$chain2[,"N"], samplesList$chain3[,"N"])), main = "Posterior N")
  col = col.list[this.mod]
  polygon(density(cbind(samplesList$chain1[,"N"], samplesList$chain2[,"N"], samplesList$chain3[,"N"])), col = col, border = col)
  abline(v = 30,lwd = 2, col = "black")
  
}

for(this.mod in 1:4){
  load(paste("~/PROJECTS/OpenRepo/NIMOUT/NimOutForSimDataSet1Rep4mod", this.mod, ".RData", sep = ""))
  if(this.mod == 1)plot(density(cbind(samplesList$chain1[,"sigma"], samplesList$chain2[,"sigma"], samplesList$chain3[,"sigma"])), main ="Posterior Sigma", xlim = c(1,10), ylim = c(0,0.5))
  if(this.mod!=1)lines(density(cbind(samplesList$chain1[,"sigma"], samplesList$chain2[,"sigma"], samplesList$chain3[,"sigma"])), main = "Posterior Sigma")
  col = col.list[this.mod]
  polygon(density(cbind(samplesList$chain1[,"sigma"], samplesList$chain2[,"sigma"], samplesList$chain3[,"sigma"])), col = col, border = col)
  abline(v = 1,lwd = 2, col = "black")
  
}


