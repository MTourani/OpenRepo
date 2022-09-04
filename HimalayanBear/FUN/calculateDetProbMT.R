#' @title NIMBLE Function to calculate the sqaure distance between an AC location and several detectors.
#'
#' @description
#' \code{calculateDistance} is a NIMBLE Function to calculate the sqaure distance between an AC location and several detectors.
#' 
#' @param sxy \code{Vector} of length n.detectors containing observation/non-observations. 
#' @param detector.xy \code{Matrix} of length n.detectors denoting the detector-specific detection probability.
#'
#' @examples
#' d2[i,1:n.detectors,t] <- calculateDistance(sxy[i,1:2,t], detector.xy[1:n.detectors,1:2,t])

calculateDetProb <- nimbleFunction(run = function( sxy = double(1)
                                                   , detector.xy = double(2)
                                                   , pZero = double(0)
                                                   , sigma = double(0)
                                                   , maxDist= double(0, default = 0.0)){
   # Return type declaration
   returnType(double(1))
   
   # Check input dimensions
   n.detectors <- dim(detector.xy)[1]
   alpha <- -1 / (2 * sigma * sigma)
   
   # Calculate distance vector (Output)
   d2 <- pow(detector.xy[1:n.detectors,1] - sxy[1], 2) + pow(detector.xy[1:n.detectors,2] - sxy[2], 2)
   
   # Calculate detector-specific probabilities
   if(maxDist > 0){
      p <- rep(0.0, n.detectors)
      for(j in 1:n.detectors){
         if(d2[j] <= (maxDist*maxDist)){p[j] <- pZero * exp(alpha * d2[j])}
         }#j
      }else{
         p <- pZero * exp(alpha * d2)
         }
   
   return(p)
   })
