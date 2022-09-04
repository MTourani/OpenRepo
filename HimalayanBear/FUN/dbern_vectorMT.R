#' @title Function to create a NIMBLE custom distribution with internalized detection
#'  probabilities calculation for faster SCR model runs.
#'
#' @description
#' \code{dSCR} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} of length n.detectors containing observation/non-observations 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param detector.xy A \code{Matrix}  of dimensions n.detectors*2 with detectors x and y coordinates.
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dSCR(p0 , sigma, sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t], z[i,t] == 2)

#### 1.Density function ####
dbern_vector <- nimbleFunction(run = function( x = double(1)
                                             , p = double(1)
                                             , indicator = double(0, default = 1.0)
                                             , log = integer(0, default = 0)){
   # Return type declaration
   returnType(double(0))
   outProb <- 0.0
   
   ## Check input dimensions
   n.detectors <- length(x)
   if(n.detectors <= 0){stop("invalid number of detectors")}

   ## Shortcut if individual is not available for detection
   if(indicator == 0){
      if(sum(x[1:n.detectors]) == 0){
         if(log == 0) return(1.0) else return(0.0)
      }else{
       if(log == 0) return(0.0) else return(-Inf)
      }
   }

   ## Calculate the likelihood of the detection observations
   detectLogLikeli <- rep(-Inf, n.detectors)
   for(j in 1:n.detectors){
      # Calculate the log-likelihood of each detection observation
      if(x[j] == 0){detectLogLikeli[j] <- log(1.0 - p[j])} else {detectLogLikeli[j] <- log(p[j])}
      # If the probability of detecting the current cell is zero then stop calculating and return the zero probability
      if(detectLogLikeli[j] <= -Inf){if(log == 0){return(0.0)}else{return(-Inf)}}
      }#j
   
   # Output
   outProb <- sum(detectLogLikeli[1:n.detectors])
   if(log == 0){outProb <- exp(outProb)}
   return(outProb)
})

#### 2.Sampling function ####
rbern_vector <- nimbleFunction(run = function( n = integer(0)
                                             , p = double(1)
                                             , indicator = double(0, default = 1.0)){
   # Return type declaration
   returnType(double(1))
   
   # Check input dimensions
   if(n!=1){print("rinhomPP only allows n = 1; using n = 1")}
   n.detectors <- length(p)
   if(n.detectors <= 0){stop("invalid number of detectors")}

   # Shortcut if individual is not available for detection
   if(indicator == 0){return(rep(0.0, n.detectors))}
 
   # Draw from a Bernoulli distribution with the calculated probability
   detectOut <- rbinom(n.detectors, 1, p)
   
   # Output
   return(detectOut)
})

#### 3.Registration ####
registerDistributions(list(
   dbern_vector = list(
      BUGSdist = "dbern_vector(p, indicator)",
      types = c( "value = double(1)", "p = double(1)", "indicator = double(0)"),
      pqAvail = FALSE)))

