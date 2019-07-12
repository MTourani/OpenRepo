
## ------ IMPORT REQUIRED LIBRARIES ------
.libPaths(c("C:/Rfolder", .libPaths())) 

library(nimble)                    
library(R.utils)                   
library(coda)

## ------ Source FUNCTIONS ------
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
                         , niter = 5000
                         , nburnin = 2000
                         , nchains = 3
                         , samplesAsCodaMCMC = TRUE)
  
  outname <- paste(path.out,"NimOutFor",set, sep = "")
  save(samplesList, file = outname)
  
  file.list <- list.files(path.in)
}

## -------------------------------------------------------------------------------------------------------------
## ----- VII.PROCESS THE OUTPUT -----  
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


