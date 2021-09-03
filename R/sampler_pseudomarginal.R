correctPercentage <- function(selectionMatrix,groundTruth){
  n= dim(selectionMatrix)[1]
  m= dim(selectionMatrix)[2]
  correctPercentage = 0
  for(i in 1:n){
    for(j in 1:m){
      if( selectionMatrix[i, j] == groundTruth[i, j]){
        correctPercentage = correctPercentage + 1
      }
    }
  }
  correctPercentage / prod(dim(groundTruth)) * 100
}

chainAutoCorrelation <- function(MHChainOfConfigurations){
  autoCorrelationMatrix = matrix(NA, 20, 20)
  for (i in 1:20) {
    for(j in 1:20){
      autoCorrelationMatrix[i, j]=acf(MHChainOfConfigurations[i, j,], lag.max = 1, plot = F)$acf[2, 1, 1]
    }
  }
  autoCorrelationMatrix
}

mcmcseEssForImage <- function(MHChainOfConfigurations){

  essMatrix = matrix(NA, 20, 20)

  for (i in 1:20) {
    for (j in 1:20){
      essMatrix[i,j]= ess(MHChainOfConfigurations[i,j,])
      cat("\n at i:",i, " and j:",j)
    }

  }

  essMatrix
}


#rename this to be petImageGenerator
#'@export
imageGeneratorPet <- function(groundTruth,
                              phi,
                              theta,
                              noiseLevel = 0.64,
                              timeFrames){
  if (length(phi) > 2) {stop("THIS FUNCTION DOES NOT WORK FOR 3 COMPARTMENTAL MODELS YET!!")}
  phi3 = c(0.0044406807, 0.0001010392,0.0014582801)
  theta3 = c(0.0004518293, 0.0027728739,  0.0107752968)
  dataMatrix = array(NA,c(n,m,32))

  for(i in 1:n){
    for(j in 1:m){
      if(groundTruth[i,j] == 1){
        dataMatrix[i,j,] =  dataSimltr(timeFrames,phi[1],theta[1],noiseLevel = noiseLevel)$Y
      }else if(groundTruth[i,j] == 2){
        dataMatrix[i,j,] = dataSimltr(timeFrames,phi,theta,noiseLevel = noiseLevel)$Y
      }else if(groundTruth[i,j] == 3){
        dataMatrix[i,j,] = dataSimltr(timeFrames,phi3,theta3,noiseLevel = noiseLevel)$Y
      }
    }
    # print(i)
  }
  cat("Image has been generated")
  dataMatrix
}

## Full Node-Wise-Pseudo-Marginal-----------------------------------------------
#'@export
nodeWisePMHSamplerPETRcpp <- function(imageData,
                                      numberOfIterations,
                                      J,
                                      pottsStateSpace = c(1,2),
                                      smcParameters = list(numSamples = 200,
                                                           numDistrbtns = 500)
){
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n",n,"m:",m)

  # currModelOrderMatrix = isingibbs2Pet(niter = 10,n,m)

  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,
                                           n = n,
                                           m = m,
                                           J=J,
                                           stateSpace = pottsStateSpace)
  print(currModelOrderMatrix)
  # stop("STAPPH")

  numStates = length(pottsStateSpace)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)


  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
  # singleMarginalEstimateMatrix = array(NA,c(n,m,2))


  cat("\n Computing marginal likelihoods using smc...")

  # cat("\n Marginal Likelihoods computed....\n")

  for(i in 1:n){
    for(j in 1:m){

      cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

      smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
      currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
      MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
                                         normalisingConstant = smcOuput$normalisingConstant,
                                         particles = smcOuput$particles[,smcParameters$numDistrbtns,],
                                         logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
    }
  }
  print(currMarginalLikelihoodMatrix)
  # print(singleMarginalEstimateMatrix[,,1])
  # print(singleMarginalEstimateMatrix[,,2])

  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        # proposalStateSpace = pottsStateSpace[!pottsStateSpace == currModelOrderMatrix[i,j]]

        ## We want to propose new states (i.e. other than the current state)
        proposalStateSpace = pottsStateSpace[!pottsStateSpace == currModelOrderMatrix[i,j]]
        ## Propose from the propossalStateSpace
        if(numStates == 2){
          proposedNodeModelOrder = proposalStateSpace[1]
        }else{
          proposedNodeModelOrder = sample(proposalStateSpace,1)
        }

        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        # proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        proposedPottsConfigPrior = pottsPriorSimplified(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates
        # currMargLik = mean(exp(unNormLogWeights[i,j,modelMatrix[i,j], ]))
        # propMargLik  = mean(exp(unNormLogWeights[i,j,propModelOrd, ] ))

        # currMargLik = dnorm(modelMatrix[i,j],0,priorSigma)*dnorm(dataMatrix[i,j],modelMatrix[i,j],likeSigma)+rnorm(1,0,1)
        #propMargLik = dnorm(propModelOrd,0,priorSigma)*dnorm(dataMatrix[i,j],propModelOrd,likeSigma)+rnorm(1,0,1)

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        smcOutput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(proposedNodeModelOrder*2)+1)
        proposedNodeMarginalLikeli = smcOutput$normalisingConstant

        cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
            "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
            # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
            "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
                                             normalisingConstant = smcOuput$normalisingConstant,
                                             particles = smcOuput$particles[,smcParameters$numDistrbtns,],
                                             logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected
          MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
      }
    }
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations, MHChainOfEstimates = MHChainOfEstimates)
}

nodeWiseMarginalModelPet <- function(MHChainOfConfigurations){

  n = dim(MHChainOfConfigurations)[1]
  m = dim(MHChainOfConfigurations)[2]
  numOfIter = dim(MHChainOfConfigurations)[3]

  positiveCounts = negativeCounts  = selectionMatrix = matrix(0,n,m)

  for(k in 2:numOfIter){
    for(i in 1:n){
      for(j in 1:m){
        # print(configArray[i,j,k])
        if(MHChainOfConfigurations[i,j,k] == 1){
          positiveCounts[i,j] = positiveCounts[i,j] + 1
        }else{
          negativeCounts[i,j]  = negativeCounts[i,j] +1
        }
      }
    }
  }

  for(i in 1:n){
    for(j in 1:m){

      if(positiveCounts[i,j] <= negativeCounts[i,j]){
        selectionMatrix[i,j] = 2
      }else{
        selectionMatrix[i,j] = 1
      }

    }
  }
  selectionMatrix
}



####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
nodeWisePMHSamplerPET <- function(imageData,
                                  numberOfIterations,
                                  J = 0.6,
                                  pottsStateSpace = 0,
                                  smcParameters = list(numSamples = 200,
                                                        numDistrbtns = 500)
){

####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n",n,"m:",m)

  currModelOrderMatrix = isingibbs2Pet(niter = 10,n,m)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  singleMarginalEstimateMatrix = array(NA,c(n,m,2))

  cat("\n Computing marginal likelihoods using smc...")

  cat("\n Marginal Likelihoods computed....\n")

  for(i in 1:n){
    for(j in 1:m){

      cat("\n Currently at:",j+(4*2*(i-1)),"/",n*m, "\n")
      if(MHChainOfConfigurations[i,j,1] == 1){
        currMarginalLikelihoodMatrix[i,j] = smcSampler(petLogLikelihoodDensity,
                                                              petLogBioPriorDensity,
                                                              petBioPriorSampler,
                                                              150,400,dimParameter = 3,
                                                              datum = imageData[i,j,],
                                                              timeFrames = timeFrames)$normalisingConstant
      }else{
        currMarginalLikelihoodMatrix[i,j] = smcSampler(petLogLikelihoodDensity,
                                                              petLogBioPriorDensity,
                                                              petBioPriorSampler,
                                                              150,400,dimParameter = 5,
                                                              datum = imageData[i,j,],
                                                              timeFrames = timeFrames)$normalisingConstant
      }

    }
  }
  print(currMarginalLikelihoodMatrix)
  print(singleMarginalEstimateMatrix[,,1])
  print(singleMarginalEstimateMatrix[,,2])

  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 1){2}else{1}
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates
        # currMargLik = mean(exp(unNormLogWeights[i,j,modelMatrix[i,j], ]))
        # propMargLik  = mean(exp(unNormLogWeights[i,j,propModelOrd, ] ))

        # currMargLik = dnorm(modelMatrix[i,j],0,priorSigma)*dnorm(dataMatrix[i,j],modelMatrix[i,j],likeSigma)+rnorm(1,0,1)
        #propMargLik = dnorm(propModelOrd,0,priorSigma)*dnorm(dataMatrix[i,j],propModelOrd,likeSigma)+rnorm(1,0,1)

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])

        proposedNodeMarginalLikeli = smcSampler(petLogLikelihoodDensity,
                                                       petLogBioPriorDensity,
                                                       petBioPriorSampler,
                                                       150,400,dimParameter = (proposedNodeModelOrder*2)+1,
                                                       datum = imageData[i,j,],
                                                       timeFrames = timeFrames)$normalisingConstant
        cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
            "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
            # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
            "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
        }else{
          #Rejected
          if(DEBUG2){cat("\n Rejected")}
        }
      }
    }
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations)
}

# Single Estimate Approximations
####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
nodeWisePMHApproxSamplerPET <- function(imageData,
                                        numberOfIterations,
                                        J = 0.6,
                                        pottsStateSpace = 0,
                                        smcParameters = list(numSamples = 200,
                                                             numDistrbtns = 500)
){
  ####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
  ####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
  ####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
  ####### Obselete---- Avoid using, nodeWisePMHSamplerPETRcpp instead -#########
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n",n,"m:",m)

  currModelOrderMatrix = isingibbs2Pet(niter = 10,n,m)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  singleMarginalEstimateMatrix = array(NA,c(n,m,2))

  cat("\n Computing marginal likelihoods using smc...")

  cat("\n Marginal Likelihoods computed....\n")

  for(i in 1:n){
    for(j in 1:m){

      for( k in 1:2){
        singleMarginalEstimateMatrix[i,j,k] = smcSampler(petLogLikelihoodDensity,
                                                         petLogBioPriorDensity,
                                                         petBioPriorSampler,
                                                         200,500,dimParameter = (k*2)+1,
                                                         datum = imageData[i,j,],
                                                         timeFrames = timeFrames)$normalisingConstant
        cat("\n Currently at:",k+(2*(j-1))+(4*2*(i-1)),"/",n*m*2, "\n")
      }

      if(MHChainOfConfigurations[i,j,1] == 1){
        currMarginalLikelihoodMatrix[i,j] = singleMarginalEstimateMatrix[i,j,1]
      }else{
        currMarginalLikelihoodMatrix[i,j] = singleMarginalEstimateMatrix[i,j,2]
      }

    }
  }
  print(currMarginalLikelihoodMatrix)
  print(singleMarginalEstimateMatrix[,,1])
  print(singleMarginalEstimateMatrix[,,2])

  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 1){2}else{1}
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates
        # currMargLik = mean(exp(unNormLogWeights[i,j,modelMatrix[i,j], ]))
        # propMargLik  = mean(exp(unNormLogWeights[i,j,propModelOrd, ] ))

        # currMargLik = dnorm(modelMatrix[i,j],0,priorSigma)*dnorm(dataMatrix[i,j],modelMatrix[i,j],likeSigma)+rnorm(1,0,1)
        #propMargLik = dnorm(propModelOrd,0,priorSigma)*dnorm(dataMatrix[i,j],propModelOrd,likeSigma)+rnorm(1,0,1)

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        if(proposedNodeModelOrder == 1){
          proposedNodeMarginalLikeli = singleMarginalEstimateMatrix[i,j,1]
        }else{
          proposedNodeMarginalLikeli = singleMarginalEstimateMatrix[i,j,2]
        }
        cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
            "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
            "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
            "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
        }else{
          #Rejected
          if(DEBUG2){cat("\n Rejected")}
        }
      }
    }
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations,
       singleMarginalEstimateMatrix = singleMarginalEstimateMatrix)
}
