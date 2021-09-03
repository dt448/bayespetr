imageGeneratorToy <- function(groundTruth){
  dataMatrix = meanMatrix = matrix(NA,n,m)

  for(i in 1:n){
      for(j in 1:m){
        meanMatrix[i,j] = rnorm(1,groundTruth[i,j],priorSigma)
        dataMatrix[i,j] = rnorm(1,meanMatrix[i,j],likeSigma)
      }
  }

  list(dataMatrix = dataMatrix, meanMatrix = meanMatrix)

}


### Full NWPM sampler for TOY using Rcpp
#'@export
nodeWisePMHSamplerToyRcpp <- function(imageData,
                                      numberOfIterations,
                                      J,
                                      pottsStateSpace,
                                      smcParameters = list(numSamples = 200,
                                                           numDistrbtns = 500),
                                      priorSigma = 5, likeSigma = 1
){
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n:", n, "m:", m)

  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,n = n,m = m,J = J,stateSpace = pottsStateSpace)

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

      for(state in 1:numStates){
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){

          smcOuput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                       numDistributions = smcParameters$numDistrbtns,
                                       datum = imageData[i,j],
                                       dimParameter = 1,
                                       priorSigma = priorSigma,
                                       likeSigma = likeSigma,
                                       mu = pottsStateSpace[state])

          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant

          MHChainOfEstimates[[i,j,1]] = list(modelOrder  = pottsStateSpace[state],
                                             normalisingConstant = smcOuput$normalisingConstant,
                                             particles = smcOuput$particles[,smcParameters$numDistrbtns],
                                             logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }
      }
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

        # proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 5){-5}else{5}
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
        smcOutput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                      numDistributions = smcParameters$numDistrbtns,
                                      datum = imageData[i,j],
                                      dimParameter = 1,
                                      priorSigma = priorSigma,
                                      likeSigma = likeSigma,
                                      mu = proposedNodeModelOrder)
        # smcOutput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
        #                        priorSampler = toyBioPriorSampler,
        #                       logPriorDensity = toyLogPriorDensity,
        #                       smcParameters$numSamples,smcParameters$numDistrbtns,
        #                       imageData[i,j],likeSigma  = likeSigma,
        #                       priorSigma = priorSigma, dimParameter = 1,mu = proposedNodeModelOrder)

        proposedNodeMarginalLikeli = smcOutput$normalisingConstant

        # cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
        # "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
                                             normalisingConstant = smcOuput$normalisingConstant,
                                             particles = smcOuput$particles[,smcParameters$numDistrbtns],
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


### NWPM-SE approximation  sampler for TOY using Rcpp
#'@export
toy_NWPM_SE <- function(imageData,
                        numberOfIterations,
                        J,
                        pottsStateSpace,
                        smcParameters = list(numSamples = 200,
                                             numDistrbtns = 500),
                        priorSigma = 5, likeSigma = 1,
                        SAVE_DIR = NULL
){
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n:", n, "m:", m)

  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,n = n,m = m,J = J,stateSpace = pottsStateSpace)

  numStates = length(pottsStateSpace)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  SE_marginal_estimates  = array(NA,c(n,m,numStates))
  SE_smc_outputs = array(list(NULL), c(n,m,numStates))

  cat("\n Computing marginal likelihoods using smc...")
  # computing the marginal estimates for all model orders for all nodes of the image
  for(state in 1:numStates){
    for(i in 1:n){
      for(j in 1:m){

        cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
        smcOuput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                                                       numDistributions = smcParameters$numDistrbtns,
                                                                       datum = imageData[i,j],
                                                                       dimParameter = 1,
                                                                       priorSigma = priorSigma,
                                                                       likeSigma = likeSigma,
                                                                       mu = pottsStateSpace[state])

        # marginal estiamtes stored in this matrix
        SE_marginal_estimates[i,j,state] = smcOuput$normalisingConstant

        # other smc outputs such as V_D stored in this array for inference
        SE_smc_outputs[[i,j,state]] = list(modelOrder  = pottsStateSpace[state],
                                               normalisingConstant = smcOuput$normalisingConstant,
                                               particles = smcOuput$particles[,smcParameters$numDistrbtns],
                                               logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])

        # store the marginal likelihood estimaets of the current state (intialised using gibbs) in this matrix.
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
        }

      }
    }
  }

  print(currMarginalLikelihoodMatrix)
  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){

        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        # proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 5){-5}else{5}
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
        # smcOutput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
        #                               numDistributions = smcParameters$numDistrbtns,
        #                               datum = imageData[i,j],
        #                               dimParameter = 1,
        #                               priorSigma = priorSigma,
        #                               likeSigma = likeSigma,
        #                               mu = proposedNodeModelOrder)
        # smcOutput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
        #                        priorSampler = toyBioPriorSampler,
        #                       logPriorDensity = toyLogPriorDensity,
        #                       smcParameters$numSamples,smcParameters$numDistrbtns,
        #                       imageData[i,j],likeSigma  = likeSigma,
        #                       priorSigma = priorSigma, dimParameter = 1,mu = proposedNodeModelOrder)

        proposedNodeMarginalLikeli = SE_marginal_estimates[i,j,which(proposedNodeModelOrder==pottsStateSpace)]

        # cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
        # "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
      }
    }
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    if(!is.null(SAVE_DIR)&& k %% 1000 == 0){
      saveRDS(MHChainOfConfigurations,paste(SAVE_DIR,"/MHChainOfConfigurations_PROGRESS.RDS",sep=""))
    }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations, SE_marginal_estimates = SE_marginal_estimates, SE_smc_outputs = SE_smc_outputs)
}


### NWPM-MA  sampler for TOY using Rcpp (BASED ON NWPM-SE above)
#'@export
toy_NWPM_MA <- function(imageData,
                        numberOfIterations,
                        J,
                        pottsStateSpace,
                        smcParameters = list(numSamples = 200,
                                             numDistrbtns = 500),
                        priorSigma = 5, likeSigma = 1,
                        SAVE_DIR = NULL,
                        REFRESH=F,refreshConst
){
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n:", n, "m:", m)

  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,n = n,m = m,J = J,stateSpace = pottsStateSpace)

  numStates = length(pottsStateSpace)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  SE_marginal_estimates  = array(NA,c(n,m,numStates))
  SE_smc_outputs = array(list(NULL), c(n,m,numStates))

  cat("\n Computing marginal likelihoods using smc...")
  # computing the marginal estimates for all model orders for all nodes of the image
  for(state in 1:numStates){
    for(i in 1:n){
      for(j in 1:m){

        cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
        smcOuput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                     numDistributions = smcParameters$numDistrbtns,
                                     datum = imageData[i,j],
                                     dimParameter = 1,
                                     priorSigma = priorSigma,
                                     likeSigma = likeSigma,
                                     mu = pottsStateSpace[state])

        # marginal estiamtes stored in this matrix
        SE_marginal_estimates[i,j,state] = smcOuput$normalisingConstant

        # other smc outputs such as V_D stored in this array for inference
        SE_smc_outputs[[i,j,state]] = list(modelOrder  = pottsStateSpace[state],
                                           normalisingConstant = smcOuput$normalisingConstant,
                                           particles = smcOuput$particles[,smcParameters$numDistrbtns],
                                           logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])

        # store the marginal likelihood estimaets of the current state (intialised using gibbs) in this matrix.
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
        }

      }
    }
  }

  print(currMarginalLikelihoodMatrix)
  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){

        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        # proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 5){-5}else{5}
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
        # smcOutput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
        #                               numDistributions = smcParameters$numDistrbtns,
        #                               datum = imageData[i,j],
        #                               dimParameter = 1,
        #                               priorSigma = priorSigma,
        #                               likeSigma = likeSigma,
        #                               mu = proposedNodeModelOrder)
        # smcOutput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
        #                        priorSampler = toyBioPriorSampler,
        #                       logPriorDensity = toyLogPriorDensity,
        #                       smcParameters$numSamples,smcParameters$numDistrbtns,
        #                       imageData[i,j],likeSigma  = likeSigma,
        #                       priorSigma = priorSigma, dimParameter = 1,mu = proposedNodeModelOrder)

        proposedNodeMarginalLikeli = SE_marginal_estimates[i,j,which(proposedNodeModelOrder==pottsStateSpace)]

        # cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
        # "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
      }
    }
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    if(!is.null(SAVE_DIR)){
      saveRDS(MHChainOfConfigurations[ , ,1:k],paste(SAVE_DIR,"MHChainOfConfigurations_PROGRESS.RDS",sep=""))
    }
    if(REFRESH == T && k %% refreshConst == 0 && k!=numberOfIterations){
      for(state in 1:numStates){
        for(i in 1:n){
          for(j in 1:m){

            cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
            smcOuput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                         numDistributions = smcParameters$numDistrbtns,
                                         datum = imageData[i,j],
                                         dimParameter = 1,
                                         priorSigma = priorSigma,
                                         likeSigma = likeSigma,
                                         mu = pottsStateSpace[state])

            # marginal estiamtes stored in this matrix
            SE_marginal_estimates[i,j,state] = smcOuput$normalisingConstant

            # other smc outputs such as V_D stored in this array for inference
            SE_smc_outputs[[i,j,state]] = list(modelOrder  = pottsStateSpace[state],
                                               normalisingConstant = smcOuput$normalisingConstant,
                                               particles = smcOuput$particles[,smcParameters$numDistrbtns],
                                               logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])

            # store the marginal likelihood estimaets of the current state (intialised using gibbs) in this matrix.
            if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
              currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
            }

          }
        }
      }
    }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations, SE_marginal_estimates = SE_marginal_estimates, SE_smc_outputs = SE_smc_outputs)
}


### NWPM-MA  sampler for TOY using Rcpp (BASED ON NWPM-SE above)
#'@export
toy_NWPM_MA2 <- function(imageData,
                        numberOfIterations,
                        J,
                        pottsStateSpace,
                        smcParameters = list(numSamples = 200,
                                             numDistrbtns = 500),
                        priorSigma = 5, likeSigma = 1,
                        SAVE_DIR = NULL,
                        REFRESH=F,refreshConst
){
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  cat("n:", n, "m:", m)

  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,n = n,m = m,J = J,stateSpace = pottsStateSpace)

  numStates = length(pottsStateSpace)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  SE_marginal_estimates  = array(NA,c(n,m,numStates))
  SE_smc_outputs = array(list(NULL), c(n,m,numStates))

  cat("\n Computing marginal likelihoods using smc...")
  # computing the marginal estimates for all model orders for all nodes of the image
  for(state in 1:numStates){
    for(i in 1:n){
      for(j in 1:m){

        cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
        smcOuput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                     numDistributions = smcParameters$numDistrbtns,
                                     datum = imageData[i,j],
                                     dimParameter = 1,
                                     priorSigma = priorSigma,
                                     likeSigma = likeSigma,
                                     mu = pottsStateSpace[state])

        # marginal estiamtes stored in this matrix
        SE_marginal_estimates[i,j,state] = smcOuput$normalisingConstant

        # other smc outputs such as V_D stored in this array for inference
        SE_smc_outputs[[i,j,state]] = list(modelOrder  = pottsStateSpace[state],
                                           normalisingConstant = smcOuput$normalisingConstant,
                                           particles = smcOuput$particles[,smcParameters$numDistrbtns],
                                           logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])

        # store the marginal likelihood estimaets of the current state (intialised using gibbs) in this matrix.
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
        }

      }
    }
  }

  print(currMarginalLikelihoodMatrix)
  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){

        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        # proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 5){-5}else{5}
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
        # smcOutput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
        #                               numDistributions = smcParameters$numDistrbtns,
        #                               datum = imageData[i,j],
        #                               dimParameter = 1,
        #                               priorSigma = priorSigma,
        #                               likeSigma = likeSigma,
        #                               mu = proposedNodeModelOrder)
        # smcOutput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
        #                        priorSampler = toyBioPriorSampler,
        #                       logPriorDensity = toyLogPriorDensity,
        #                       smcParameters$numSamples,smcParameters$numDistrbtns,
        #                       imageData[i,j],likeSigma  = likeSigma,
        #                       priorSigma = priorSigma, dimParameter = 1,mu = proposedNodeModelOrder)

        proposedNodeMarginalLikeli = SE_marginal_estimates[i,j,which(proposedNodeModelOrder==pottsStateSpace)]

        # cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
        # "currPottsConfigPrior: ",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
        # cat("\n acceptance Ration:",accRatio)
        if(accRatio >= runif(1)){
          if(DEBUG2){cat("\n Accepted")}
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
      }
    }
    #record whole image for that iteration
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    if(REFRESH == T && k %% refreshConst == 0 && k!=numberOfIterations){
      MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
      if(!is.null(SAVE_DIR)){
        saveRDS(MHChainOfConfigurations,paste(SAVE_DIR,"/MHChainOfConfigurations_PROGRESS.RDS",sep=""))
      }
      for(state in 1:numStates){
        for(i in 1:n){
          for(j in 1:m){

            cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
            smcOuput = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
                                         numDistributions = smcParameters$numDistrbtns,
                                         datum = imageData[i,j],
                                         dimParameter = 1,
                                         priorSigma = priorSigma,
                                         likeSigma = likeSigma,
                                         mu = pottsStateSpace[state])

            # marginal estiamtes stored in this matrix
            acceration = smcOuput$normalisingConstant/SE_marginal_estimates[i,j,state]
            if( acceration >= runif(1)){
              SE_marginal_estimates[i,j,state] = smcOuput$normalisingConstant

              # other smc outputs such as V_D stored in this array for inference
              SE_smc_outputs[[i,j,state]] = list(modelOrder  = pottsStateSpace[state],
                                                 normalisingConstant = smcOuput$normalisingConstant,
                                                 particles = smcOuput$particles[,smcParameters$numDistrbtns],
                                                 logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])

              # store the marginal likelihood estimaets of the current state (intialised using gibbs) in this matrix.
              if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
                currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
              }
            }
          }
        }
      }
    }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations, SE_marginal_estimates = SE_marginal_estimates, SE_smc_outputs = SE_smc_outputs)
}





nodeWisePMHApproxSampler <- function(imageData,
                                     numberOfIterations,
                                     J = 0.5,
                                     pottsStateSpace = 0,
                                     smcParameters = list(numSamples = 200, numDistrbtns=15),
                                     toyModelParameters = list(likeSigma = 1, priorSigma = 5)){
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  currModelOrderMatrix = isingibbs2(niter = 10,n,m)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  singleMarginalEstimateMatrix = array(NA,c(n,m,2))

  for(i in 1:n){
    for(j in 1:m){
      for(k in 1:2){
        if(k == 1){
          singleMarginalEstimateMatrix[i,j,k] = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                                                           logPriorDensity = toyLogPriorDensity,
                                                           priorSampler = toyBioPriorSampler,
                                                           numSamples = smcParameters$numSamples,
                                                           numDistributions = smcParameters$numDistrbtns,
                                                           dimParameter = 1,
                                                           datum = imageData[i,j],
                                                           priorSigma = toyModelParameters$priorSigma,
                                                           likeSigma = toyModelParameters$likeSigma,
                                                           mu = 5)$normalisingConstant
        }else{
          singleMarginalEstimateMatrix[i,j,k] = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                                                           logPriorDensity = toyLogPriorDensity,
                                                           priorSampler = toyBioPriorSampler,
                                                           numSamples = smcParameters$numSamples,
                                                           numDistributions = smcParameters$numDistrbtns,
                                                           dimParameter = 1,
                                                           datum = imageData[i,j],
                                                           priorSigma = toyModelParameters$priorSigma,
                                                           likeSigma = toyModelParameters$likeSigma,
                                                           mu = -5)$normalisingConstant
        }
      }
    }
  }


  for(i in 1:n){
    for(j in 1:m){
      if(MHChainOfConfigurations[i,j,1] == 5){
        currMarginalLikelihoodMatrix[i,j] = singleMarginalEstimateMatrix[i,j,1]
      }else{
        currMarginalLikelihoodMatrix[i,j] = singleMarginalEstimateMatrix[i,j,2]
      }
    }
  }

  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 5){-5}else{5}
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates
        # currMargLik = mean(exp(unNormLogWeights[i,j,modelMatrix[i,j], ]))
        # propMargLik  = mean(exp(unNormLogWeights[i,j,propModelOrd, ] ))

        # currMargLik = dnorm(modelMatrix[i,j],0,priorSigma)*dnorm(dataMatrix[i,j],modelMatrix[i,j],likeSigma)+rnorm(1,0,1)
        #propMargLik = dnorm(propModelOrd,0,priorSigma)*dnorm(dataMatrix[i,j],propModelOrd,likeSigma)+rnorm(1,0,1)

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        if(proposedNodeModelOrder == 5){
          proposedNodeMarginalLikeli = singleMarginalEstimateMatrix[i,j,1]
        }else{
          proposedNodeMarginalLikeli = singleMarginalEstimateMatrix[i,j,2]
        }

        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior)/(currMarginalLikelihoodMatrix[i,j]*currPottsConfigPrior)
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

nodeWiseMarginalModelToy <- function(MHChainOfConfigurations){

  n = dim(MHChainOfConfigurations)[1]
  m = dim(MHChainOfConfigurations)[2]
  numOfIter = dim(MHChainOfConfigurations)[3]

  positiveCounts = negativeCounts  = selectionMatrix = matrix(0,n,m)

  for(k in 2:numOfIter){
    for(i in 1:n){
      for(j in 1:m){
        # print(configArray[i,j,k])
        if(MHChainOfConfigurations[i,j,k] == 5){
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
        selectionMatrix[i,j] = -5
      }else{
        selectionMatrix[i,j] = 5
      }

    }
  }
  selectionMatrix
}


##### Functions for Testing####################################################
#######OBSELETE DO NOT USE ####################################################
nodeWisePMHSamplerToy <- function(imageData,
                             numberOfIterations,
                             J = 0.5,
                             pottsStateSpace = 0,
                             smcParameters = list(numSamples = 100, numDistrbtns=15),
                             toyModelParameters = list(likeSigma = 1, priorSigma = 5)){
#######OBSELETE DO NOT USE ####################################################
#######OBSELETE DO NOT USE ####################################################
#######OBSELETE DO NOT USE ####################################################
#######OBSELETE DO NOT USE ####################################################
#######OBSELETE DO NOT USE ####################################################
  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]

  currModelOrderMatrix = isingibbs2(niter = 10,n,m)
  currMarginalLikelihoodMatrix = matrix(NA,n,m)

  MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
  MHChainOfConfigurations[,,1] = currModelOrderMatrix

  for(i in 1:n){
    for(j in 1:m){
       if(MHChainOfConfigurations[i,j,1] == 5){
        currMarginalLikelihoodMatrix[i,j] = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                          logPriorDensity = toyLogPriorDensity,
                          priorSampler = toyBioPriorSampler,
                          numSamples = smcParameters$numSamples,
                          numDistributions = smcParameters$numDistrbtns,
                          dimParameter = 1,
                          datum = imageData[i,j],
                          priorSigma = toyModelParameters$priorSigma,
                          likeSigma = toyModelParameters$likeSigma,
                          mu = 5)$normalisingConstant
      }else{

        currMarginalLikelihoodMatrix[i,j] = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                          logPriorDensity = toyLogPriorDensity,
                          priorSampler = toyBioPriorSampler,
                          numSamples = smcParameters$numSamples,
                          numDistributions = smcParameters$numDistrbtns,
                          dimParameter = 1,
                          datum = imageData[i,j],
                          priorSigma = toyModelParameters$priorSigma,
                          likeSigma = toyModelParameters$likeSigma,
                          mu = -5)$normalisingConstant
      }
    }
  }

  for(k in 2:numberOfIterations){
    cat("\n Currently at:",k,"/",numberOfIterations)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 5){-5}else{5}
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
  # -------------------------------------------------------------------------



        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates
        # currMargLik = mean(exp(unNormLogWeights[i,j,modelMatrix[i,j], ]))
        # propMargLik  = mean(exp(unNormLogWeights[i,j,propModelOrd, ] ))

        # currMargLik = dnorm(modelMatrix[i,j],0,priorSigma)*dnorm(dataMatrix[i,j],modelMatrix[i,j],likeSigma)+rnorm(1,0,1)
        #propMargLik = dnorm(propModelOrd,0,priorSigma)*dnorm(dataMatrix[i,j],propModelOrd,likeSigma)+rnorm(1,0,1)

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        proposedNodeMarginalLikeli = currMarginalLikelihoodMatrix[i,j] = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                          logPriorDensity = toyLogPriorDensity,
                          priorSampler = toyBioPriorSampler,
                          numSamples = smcParameters$numSamples,
                          numDistributions = smcParameters$numDistrbtns,
                          dimParameter = 1,
                          datum = imageData[i,j],
                          priorSigma = toyModelParameters$priorSigma,
                          likeSigma = toyModelParameters$likeSigma,
                          mu = proposedNodeModelOrder)$normalisingConstant
        #Compute the acceptance ratio for MH correction.
        # cat("\n proposedNodeMarginalLikeli Ration:",proposedNodeMarginalLikeli)
        # cat("\n proposedPottsConfigPrior Ration:",proposedPottsConfigPrior)
        # cat("\n currMarginalLikelihoodMatrix Ration:",currMarginalLikelihoodMatrix)
        # cat("\n currPottsConfigPrior Ration:",currPottsConfigPrior)
        accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior)/(currMarginalLikelihoodMatrix[i,j]*currPottsConfigPrior)
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
  MHChainOfConfigurations
}
