# These functions do not need all the parameters (but due to variadastic passing,
# needs to have all the parameters)
toyBioPriorSampler <- function(dimParameter, numSamples, priorSigma, mu, likeSigma){
  rnorm(numSamples, mu, priorSigma)
}

toyLogLikelihoodDensity <- function(datum, parameters, dimParameters,priorSigma, mu, likeSigma){
 dnorm(datum, parameters, likeSigma, log = T)
}

toyLogPriorDensity <- function(parameters, dimParameters, priorSigma, mu, likeSigma){
  dnorm(parameters, mu, priorSigma, log = T)
}

toyModelTrueNormalizingConst <- function( dataset, modelOrder, priorSigma, likelihoodSigma){

  sigmaFG = (priorSigma^2 * likelihoodSigma^2) / (priorSigma^2+likelihoodSigma^2)

  scaleFactor = exp(-0.5 * ((modelOrder-dataset)^2 / ((priorSigma^2 + likelihoodSigma^2))))

  # scaleFactor

  (sqrt(2 * pi * sigmaFG) * scaleFactor)/(sqrt(2 * pi * likelihoodSigma^2) * sqrt(2 * pi * priorSigma^2))
}
