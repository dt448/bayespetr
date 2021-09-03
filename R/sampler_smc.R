mcmcMovesCESSgeneral = function(logLikelihoodDensity,
                                logPriorDensity,
                                initialValue,
                                steps = 4,
                                dimParameter,
                                alphaT,
                                datum,
                                cholDecom,
                                ...){
  # dimParameter = (modelOrder*2)+1
  modelOrder = (dimParameter-1)/2
  # if(DEBUG){cat("\n dim:",dimParameter,"\n")}

  # X = array(0,c(1,steps,dimParameters))
  X = matrix(NA,(steps+1),dimParameter)
  X [1,] = initialValue
  Y = rep(NA,dimParameter)

  accRate = 0

  for(i in 2:(steps+1)){
    propsedNoise = rnorm(dimParameter)
    # Y = X[i-1,] + t(cholDecom)%*%rnorm(dimParameter) #propsedNoise
    Y = X[i-1,] + t(cholDecom)%*%propsedNoise #propsedNoise
    # cat("\n proposedNoise:", propsedNoise, "\n Y:",Y)
    # if(DEBUG && Y[dimParameter]>=-100){ cat("\n", i ," mcmcmoves Y:",Y, "X:",X[1,i-1,],"\n")}
    # if(DEBUG2){cat("\n b:",Y)}

    propVPriorCheck<-TRUE

    propLogPriorDensity = logPriorDensity(Y,dimParameter, ...)
    if(exp(propLogPriorDensity) == 0){
      propVPriorCheck <- FALSE}
    # for(m in 1:modelOrder){
    #   if(phivPrior(Y[m]) == 0 | thetavPrior(Y[modelOrder+m]) == 0){
    #     if(DEBUG){print('\n  NEGATIVE PROP VALUE \n')}
    #     # in this case the propValue is NOT in the prior so fails the propVPriorCheck
    #     propVPriorCheck<-FALSE
    #   }
    # }

    if(propVPriorCheck == FALSE){ a = 0 }else{

      # YLogLikelihood = petLogLikelihoodDensity(datum, Y, dimParameter, timeFrames)
      # XLogLikelihood = petLogLikelihoodDensity(datum, X[i-1,], dimParameter, timeFrames)
      # currLogPriorDensity = petLogBioPriorDensity(X[i-1,],dimParameter)

      YLogLikelihood = logLikelihoodDensity(datum, Y, dimParameter, ...)
      XLogLikelihood = logLikelihoodDensity(datum, X[i-1,], dimParameter, ...)
      currLogPriorDensity = logPriorDensity(X[i-1,], dimParameter, ...)

      YLogProb = (YLogLikelihood* alphaT) + propLogPriorDensity

      XLogProb = (XLogLikelihood* alphaT) + currLogPriorDensity


      #print(YLogProb)
      #print(XLogProb)
      if(YLogProb == -Inf || is.na(YLogProb)){
        a = 0} # if propValue
      else if (XLogProb == -Inf|| is.na(YLogProb)){
        a = 1}
      else{
        a = exp(YLogProb-XLogProb )}


      # if(DEBUG){cat("\n mcmcmoves a:",a, "Yprob:",YLogProb,"Xprob:",XLogProb,"\n")}

    }


    u = runif(1)
    if(u<a){
      X[i,] <- Y
      # cat("\n Accepted")
      # if(DEBUG){cat("\n mcmc : Accepted ")}
      accRate = accRate +1
    }else{
      X[i,] <- X[i-1,]
      # cat("\n Rejected")
      # if(DEBUG){cat("\n mcmc: Rejected ")}
    }

  }
  accRate = accRate/steps
  list(X=X[(steps+1),], accRate = accRate)

}

#'SMC Sampler for Bayesian Simulated Annealing
#'
#'\code{smcSamplerGeneral} Returns a List containing various values of the SMC
#'sampler e.g.  particle and (log) weight collection, normalising constant of
#'posterior denisty and other adaptive techniques metadata.
#'
#'This is a general function that will provide weighted samples that can be used
#'to approximate different properties of the Bayesian model of interest
#'
#'@param  logLikelihoodDensity A function that returns the log denisty of the
#'  likelihood.
#'@param logPriorDensity A function that retunrs the log denisty of the prior.
#'@param priorSampler A functoin that returns samples from the prior.
#'@param numSamples A number specifiying the number of particles to be used in
#'  the SMC sampler.
#'@param numDistributions A number specifiying the number of intermediate
#'  distributions to be used in the simulated annealing (when adpatice CESS is
#'  not used).
#'@param dimParameter The dimension of the paramaters space.
#'@param datum A vector containing the data point.
#'@param ... Other varaibles to be passed onto likelihood, prior etc.
#'
#'@examples
#'
#'\dontrun{
#' smcOutput = smcSampler(logLikelihoodDensity = petLogLikelihoodDensity,
#' logPriorDensity = petLogBioPriorDensity,
#' priorSampler = petBioPriorSampler,
#' dimParameter = 5,
#' datum = data1,
#' timeFrames = timeFrames)
#'}
#'
#'@export
smcSampler = function(logLikelihoodDensity,
                      logPriorDensity,
                      priorSampler,
                      numSamples = 200,
                      numDistributions = 500,
                      dimParameter,
                      datum,
                      ...){

  if(!exists("DEBUG")){DEBUG = F}
  if(!exists("DEBUG2")){DEBUG2 = F}

  if(dimParameter == 1){
    weightedAdaptiveStepSize <- adaptiveWeightedStepSizeToy
    adaptiveStepSize <- adaptiveStepSizeToy

  }
  #------ Setting up arrays to hold particles, weights etc----------------------
  particles = array(NA,c(numSamples,numDistributions,dimParameter))
  logWeights = matrix(NA,nrow = numSamples,ncol = numDistributions)
  unnormLogWeights = incremWeights = barLogWeights = acceptanceRates = logWeights
  reSamplingIndex = matrix(NA,numSamples,1)

  # Approximate ESS, for adaptive resampling
  ESS = rep(NA,numDistributions)

  # Scaled sample covariance, for adaptive MCMC kernel.
  optimalStepSizes = array(NA,c( dimParameter,numDistributions,dimParameter) )

  #------ Seting up variables for adaptive CESS --------------------------------
  # Setting threshold for cessStar.
  cessStar = 0.66*numSamples
  if(DEBUG2){print(cessStar)}
  #condESS = rep(NA,numDistrbtns-1)
  alphaVec = rep(NA,numDistributions)
  alphaVec[1] = 0

  #------- Initialising particles and weights, using the prior sampler----------

  # Using priorSampler here to generate the first set of samples
  particles[,1,] = priorSampler(dimParameter, numSamples,...)
  # particles[,1,] = firstSample
  # Assuming that they are equally distributed
  incremWeights[,1] = log(1)
  logWeights[,1] = -log(numSamples) # initialise giving each equal weight i.e. log(1/n)

  # Computing the ESS and resampling if neccesary (Note: probably don't need to
  # to do this).

  ESS[1]  = (1/sum(exp(logWeights[,1])^2))

  # Resampling if below threshold.
  if(ESS[1]<(0.5*numSamples)){
    # cat("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Resampling \n")
    particles[,1,] = particles[sample(seq(1:numSamples),numSamples,replace = TRUE,prob = exp(logWeights[,1])),1,]
    barLogWeights[,1] = - log(numSamples)
  }else{
    barLogWeights[,1] = logWeights[,1]
  }

  #------ Begin propogation and resampling -------------------------------------
  for(d in 2:numDistributions){

    if(DEBUG2){cat("\n ############################################################################################ Currently at",d,"/",numDistributions, "\n")}

    if(alphaVec[d-1] == 1){
      break
    }

    # binSearOutput = binsearch(CESS,range = c(alphaVec[d-1],1),target = cessStar, threshold = 10,alphaTm = alphaVec[d-1],prevNormWeights = logWeights[,d-1],prevSamples = samples[,d-1,] ,y = y,timeFrames = timeFrames,modelId = modelId,numSamples=numSamples,numDistrbtns=numDistrbtns)$where
    binSearOutput = (d/numDistributions)^5

    if(length(binSearOutput) == 2){
      alphaStar =binSearOutput[2]
    }else{
      alphaStar = binSearOutput
    }

    # if(DEBUG2){cat("\n alphaStar: ",alphaStar)}

    if(alphaStar<=1){
      alphaVec[d] = alphaStar
      # if(DEBUG2){cat(", alphaT",alphaVec[d],"and alpha_{T-1}:",alphaVec[d-1])}
    }else{
      alphaVec[d] = 1
    }

    resampCheck = F

    # Compute ESS and resampling if needed.
    ESS[d] = 1/sum(exp(logWeights[,d-1])^2)
    # if(DEBUG2){cat("\n ESS:",ESS[d])}

    if( is.na(ESS[d]) ){
      # cat("\n Resampling due to NA \n")
      reeeSampleindex = sample(seq(1:numSamples),numSamples,replace = TRUE,prob = exp(logWeights[,d-1]))
      reSamplingIndex = cbind(reSamplingIndex,reeeSampleindex)
      particles[,d-1,] = particles[reeeSampleindex,d-1,]
      barLogWeights[,d-1] = -log(numSamples)
      stepSize = adaptiveStepSize(particles[,d-1,])$scaling
      resampCheck = T
    }else if( ESS[d] <(0.5*numSamples) ){
      # cat("\n Resampling due to low ESS \n")
      # particles[,d-1,] = particles[sample(seq(1:numSamples),numSamples,replace = TRUE,prob = exp(logWeights[,d-1])),d-1,]
      reeeSampleindex = sample(seq(1:numSamples),numSamples,replace = TRUE,prob = exp(logWeights[,d-1]))
      reSamplingIndex = cbind(reSamplingIndex,reeeSampleindex)
      particles[,d-1,] = particles[reeeSampleindex,d-1,]
      barLogWeights[,d-1] = -log(numSamples)
      stepSize = adaptiveStepSize(particles[,d-1,])$scaling
      # stepSize = weightedAdaptiveStepSize(samples[,d-1,],logWeights[,d-1])$scaling #Computing the optimal step size.
      resampCheck = T
    }else{
      # cat("\n not Resampling \n")
      barLogWeights[,d-1] = logWeights[,d-1]
      stepSize = weightedAdaptiveStepSize(particles[,d-1,],logWeights[,d-1])$scaling #Computing the optimal step size.
    }

    optimalStepSizes[,d,] = stepSize
    cholDecom = chol(stepSize)
    # cat("\n chol:")
    # print(cholDecom)

    for(i in 1:numSamples){

      prevSamplesLogLikeli = logLikelihoodDensity(datum,particles[i, d-1, ], dimParameter, ...)

      incremWeights[i,d] = (alphaVec[d]*prevSamplesLogLikeli) -(alphaVec[d-1]*prevSamplesLogLikeli)

      logWeights[i,d] = barLogWeights[i,d-1]+incremWeights[i,d]


      # Moving with M-H RW Kernel

      # propgParticles = mcmcMovesCESS(particles[i,d-1,],steps = 4,(dimParameter-1)/2,alphaVec[d],datum,cholDecom)

      propgParticles = mcmcMovesCESSgeneral(logLikelihoodDensity = logLikelihoodDensity,logPriorDensity = logPriorDensity,
                                            initialValue = particles[i,d-1,],steps = 3,dimParameter = dimParameter,
                                            alphaT = alphaVec[d],datum = datum,cholDecom = cholDecom, ...)
      # mcmcMovesCESSgeneral(logLikelihoodDensity,logPriorDensity,particles[i,d-1,],dimParameter = dimParameter,)
      particles[i,d,] = propgParticles$X
      # cat("\n prop:", propgParticles$X)
      # if(i==123 && d == 12){stop()}
      acceptanceRates[i,d] = propgParticles$accRate

    }


    ###################################################### Computing the normalized log Weights ################
    unnormLogWeights[,d] = logWeights[,d]
    maxLW = max(logWeights[,d])
    # if(DEBUG2){ cat("MAX WEIGHTED:",particles[match(maxLW,logWeights[,d]),d,] )}
    logWeights[,d] = logWeights[,d] - maxLW - log(sum(exp(logWeights[,d]-maxLW)))
    ############################################################################################################

  }

  prod1=mean(exp(incremWeights[,1]))
  prod1
  for(i in 2:(which(alphaVec == 1))){
    prod1 = prod1 * sum(exp(incremWeights[,i]+barLogWeights[,i-1]))

  }

  list(particles = particles, normLogWeights = logWeights, UnNormlogWeights = unnormLogWeights,
       alphaVec = alphaVec, acceptanceRates = acceptanceRates, incremWeights = incremWeights,
       barLogWeights = barLogWeights, ESS = ESS, optimalStepSizes = optimalStepSizes,
       normalisingConstant = prod1, reSamplingIndex = reSamplingIndex)
}

adaptiveStepSize = function(samples){
  dimens  = dim(samples)[2]
  print(dimens)
  sigmaaa = cov(samples)
  outPut = list(scaling = ((2.3^2)/dimens)*sigmaaa)
  #print(outPut)
  return(outPut)

}

weightedAdaptiveStepSize = function(samples,normLogWeights){
  dimens  = dim(samples)[2]
  print(dimens)
  #samples[,dimens] = log(samples[,dimens])
  #sigmaaa = cov(samples)
  sigmaaa = cov.wt(samples,exp(normLogWeights))$cov
  # print(sigmaaa)
  outPut = list(scaling = ((2.3^2)/dimens)*sigmaaa)
  #print(outPut)
  return(outPut)

}


######-- Helper function for Rcpp Implementation-------------------------------
### NOTE: This only works for PET model currently.


#'SMC Sampler for (PET Model) - Rcpp Implementation
#'
#'
#'IMPORTANT: This function only works for the PET model!
#'Returns a List containing various values of the SMC
#'sampler e.g.  particle and (log) weight collection, normalising constant of
#'posterior denisty and other adaptive techniques metadata.
#'
#'This is NOT a general function. It will return provide weighted samples that
#'can be used to approximate different properties of the PET Bayesian model.
#'
#'@param numSamples A number specifiying the number of particles to be used in
#'  the SMC sampler.
#'@param numDistributions A number specifiying the number of intermediate
#'  distributions to be used in the simulated annealing (when adpatice CESS is
#'  not used).
#'@param dimParameter The dimension of the paramaters space.
#'@param datum A vector containing the data point.
#'
#'@examples
#'
#'\dontrun{
#' smcOutput = smcSamplerRcpp( dimParameter = 5, datum = data1,
#' timeFrames = timeFrames)
#'}
#'
#'@export
smcSamplerRcpp <-function(numSamples = 200,
                          numDistributions = 500,
                          datum,
                          timeFrames,
                          dimParameter) {
  # This is to allow for unit test ( so that the Rcpp implementation is identical
  # to the R implemenatation).
  firstSample = petBioPriorSampler(dimParameter, numSamples, timeFrames)
  rTakeit(firstSample)
  smcSamplerCpp(numSamples, numDistributions, datum, timeFrames, dimParameter)

}
######--Toy Mode Worker Functions--##############################################

#'@export
smcSamplerToyRcpp <-function(numSamples = 200,
                             numDistributions = 500,
                             datum,
                             dimParameter,
                             priorSigma,
                             likeSigma,
                             mu) {
  # This is to allow for unit test ( so that the Rcpp implementation is identical
  # to the R implemenatation).
  firstSample = toyBioPriorSampler(dimParameter, numSamples, priorSigma, mu, likeSigma)
  # firstSample = toyBioPriorSampler(dimParameter, numSamples, mu, likeSigma)
  rTakeitToy(firstSample)
  smcSamplerToyCpp(numSamples, numDistributions, datum, dimParameter,priorSigma,
                   likeSigma,
                   mu)
}

adaptiveStepSizeToy = function(samples,modelId = 1){

  if(modelId != 1){sigmaaa = cov(samples)}
  else{sigmaaa = var(samples)}
  outPut = list(scaling = ((2.3^2)/modelId)*sigmaaa)
  return(outPut)

}

adaptiveWeightedStepSizeToy = function(samples,logWeights,modelId  = 1){

  # if(modelId != 1){sigmaaa = cov(samples)}
  # else{sigmaaa = var(samples)}
  wt = exp(logWeights)
  sigmaaa = sum(wt * (samples - weighted.mean(samples,wt))^2)
  outPut = list(scaling = ((2.3^2)/modelId)*sigmaaa)
  return(outPut)

}

##------------------------ Analysis

# TODO: can just use the model1V_D inputs to do the analysis
smc_indep_pet_MO <-function(model1Marginals,model2Marginals,model3Marginals,
                            model1V_D,
                            model2V_D,
                            model3V_D){
  cat("\n # TODO: can just use the model1V_D inputs to do the analysis")
  image_dim = dim(model1Marginals)
  m=image_dim[1]
  n=image_dim[2]


  selectedModelOrder = matrix(NA,m,n)

  volumeOfDistribution = matrix(NA,m,n)

  for(i in 1:m){
    for(j in 1:n){
      if(model1Marginals[i,j]>model2Marginals[i,j] && model1Marginals[i,j]>model3Marginals[i,j]){
        selectedModelOrder[i,j] = 1
        volumeOfDistribution[i,j] = model1V_D[[i,j]]$volumeOfDistribution
      }else if(model2Marginals[i,j]>model3Marginals[i,j]){
        selectedModelOrder[i,j] = 2
        volumeOfDistribution[i,j] = model2V_D[[i,j]]$volumeOfDistribution
      }else{
        selectedModelOrder[i,j] = 3
        volumeOfDistribution[i,j] = model3V_D[[i,j]]$volumeOfDistribution
      }
    }
  }
   list(selectedModelOrder = selectedModelOrder, volumeOfDistribution = volumeOfDistribution)
}


constRatesToParam = function(kVec){
  dimRateConstants = length(kVec)

  if(dimRateConstants == 2){

    phi1= kVec[1]
    theta1 = kVec[2]

    param = list(phi = phi1, theta = theta1)

  }else if(dimRateConstants == 4){

    delta = sqrt((kVec[2]+kVec[3]+kVec[4])^2-4*kVec[2]*kVec[4])

    theta1 = (kVec[2]+kVec[3]+kVec[4]+delta)/(2)
    theta2 = (kVec[2]+kVec[3]+kVec[4]-delta)/(2)

    phi1 = (kVec[1]*(theta1-kVec[3]-kVec[4]))/(delta)
    phi2 = (kVec[1]*(theta2-kVec[3]-kVec[4]))/(-delta)

    param = list(phi = c(phi1,phi2), theta = c(theta1,theta2))

  }else if(dimRateConstants == 6){

    gamma1 = kVec[2]+kVec[3]+kVec[4]+kVec[5]+kVec[6]
    gamma2 = (kVec[2]*kVec[4])+(kVec[2]*kVec[6])+(kVec[3]*kVec[6])+(kVec[4]*kVec[5])+(kVec[4]*kVec[6])
    gamma3 = kVec[2]*kVec[4]*kVec[6]

    delta1 = -(1/9)*(3*gamma2-gamma1^2)
    delta2 = (1/54)*(2*gamma1^3 - 9*gamma1*gamma2 + 27*gamma3)


    sqrtRatio = sqrt((delta2^2)/(delta1^3))


    if(delta2<0){
      upsilon = acos(sqrtRatio)
    }else if(delta2>0){
      upsilon = acos(-sqrtRatio)
    }

    theta1 = (gamma1/3) - 2*sqrt(delta1)*cos(upsilon/3)
    theta2 = (gamma1/3) - 2*sqrt(delta1)*cos((upsilon+2*pi)/3)
    theta3 = (gamma1/3) - 2*sqrt(delta1)*cos((upsilon+4*pi)/3)

    phi1 = (kVec[1]*(kVec[3]*(kVec[6]-theta1)+(kVec[4]-theta1)*(kVec[5]+kVec[6]-theta1)))/((theta1-theta2)*(theta1-theta3))
    phi2 = (kVec[1]*(kVec[3]*(kVec[6]-theta2)+(kVec[4]-theta2)*(kVec[5]+kVec[6]-theta2)))/((theta2-theta1)*(theta2-theta3))
    phi3 = (kVec[1]*(kVec[3]*(kVec[6]-theta3)+(kVec[4]-theta3)*(kVec[5]+kVec[6]-theta3)))/((theta3-theta1)*(theta3-theta2))

    param = list(phi = c(phi1,phi2,phi3), theta = c(theta1,theta2,theta3))

  }

  param
}

acos(-0.9)

volumeOfDist = function(x,rateConstants=T){
  V_D = 0
  if(rateConstants){
    firstRatio = 0
    for(i in 1:(length(x)/2)){

      if(i==1){
        firstRatio = x[1]/x[2]
        V_D = V_D + firstRatio

      }else{
        V_D = V_D + firstRatio*(x[2*i-1]/x[2*i])
      }

    }
  }else{
    for(i in 1:(length(x)/2)){
      V_D = V_D + (x[i]/x[i+(length(x)/2)])
    }
  }
  V_D
}
