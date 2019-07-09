
## ------ IMPORT REQUIRED LIBRARIES ------
.libPaths(c("C:/Rfolder", .libPaths())) 

library(nimble)                    
library(R.utils)                   
library(coda)

## ------ Source FUNCTIONS ------
# sourceDirectory("~/PROJECTS/Rovquant/Source",modifiedOnly=FALSE) 
sourceDirectory("~/PROJECTS/OpenRepo/FUN/",modifiedOnly=FALSE) 

## ------ DEFINED INPUT DIRECTORY ------
path.in <- "~/PROJECTS/OpenRepo/NIMIN/"
path.proc <- "~/PROJECTS/OpenRepo/NIMPROC/"
path.out <- "~/PROJECTS/OpenRepo/NIMOUT/"

file.list <- list.files(path.in)

## LOOP OVER INPUT FILES
while(length(file.list)>0){
  
  # Sample one file at random
  set <- sample(file.list, 1)
  
  # Load the sampled data
  load(paste(path.in, set, sep=""))
  print(set)
  file.rename(from = paste(path.in, set, sep = ""), to = paste(path.proc, set, sep = "")) 
  
  nimble.model <- nimbleModel(code = mod
                              , constants = constants
                              , data = data
                              , inits = inits
                              , check = FALSE)
  
  ### ==== 4.MODEL COMPILATION  ====
  Sys.setenv(PATH = paste("C:/Rfolder/Rtools/bin", Sys.getenv("PATH"), sep=";"))
  Sys.setenv(BINPREF = "C:/Rfolder/Rtools/mingw_$(WIN)/bin/")
  
  Cmodel <- compileNimble(nimble.model)
  mcmcSCR <- configureMCMC(nimble.model, monitors = params.mod)
  SCRMCMC <- buildMCMC(mcmcSCR)
  CompSCRMCMC <- compileNimble(SCRMCMC, project = nimble.model)
  
  samplesList <- runMCMC(CompSCRMCMC
                         , niter = 2000
                         , nburnin = 500
                         , nchains = 3
                         , samplesAsCodaMCMC = TRUE)
  
  outname <- paste(path.out,"NimOutFor",set, sep = "")
  save(samplesList, file = outname)
  
  file.list <- list.files(path.in)
}

## -------------------------------------------------------------------------------------------------------------
## ----- VII.PROCESS THE OUTPUT -----  
load("~/PROJECTS/OpenRepo/NIMOUT/NimOutForSimDataSet1Rep4mod1.RData")
summary(mcmc.list(samplesList))

for(this.mod in 1:4){
  load(paste("~/PROJECTS/OpenRepo/NIMOUT/NimOutForSimDataSet1Rep4mod",this.mod,".RData",sep=""))
  if(this.mod==1)plot(density(cbind(samplesList$chain1[,"N"],samplesList$chain2[,"N"],samplesList$chain3[,"N"])), main="Posterior N"    ,xlim=c(10,100),ylim=c(0,0.1))
  if(this.mod!=1)lines(density(cbind(samplesList$chain1[,"N"],samplesList$chain2[,"N"],samplesList$chain3[,"N"])), main="Posterior N")
  col= col.list[this.mod]
  polygon(density(cbind(samplesList$chain1[,"N"],samplesList$chain2[,"N"],samplesList$chain3[,"N"])), col=col, border=col)
  abline(v=48,lwd=2,col="black")
  
}

