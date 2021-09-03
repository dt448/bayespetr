test_that("toy NWPM works", {

  nodeWisePMHSamplerToyRcppOLD <- function(imageData,
                                           numberOfIterations,
                                           J = 0.6,
                                           pottsStateSpace = 0,
                                           smcParameters = list(numSamples = 200,
                                                                numDistrbtns = 500),
                                           priorsigma = 5, likesigma = 1
  ){
    cat("\n smcParameters:", smcParameters$numSamples)
    n = dim(imageData)[1]
    m = dim(imageData)[2]

    # cat("n",n,"m:",m)

    currModelOrderMatrix = isingibbs2Toy(niter = 10,n,m,J)
    currMarginalLikelihoodMatrix = matrix(NA,n,m)

    MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
    MHChainOfConfigurations[,,1] = currModelOrderMatrix

    MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
    # singleMarginalEstimateMatrix = array(NA,c(n,m,2))


    # cat("\n Computing marginal likelihoods using smc...")

    # cat("\n Marginal Likelihoods computed....\n")

    for(i in 1:n){
      for(j in 1:m){

        # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
        if(MHChainOfConfigurations[i,j,1] == 5){

          smcOuput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                                priorSampler = toyBioPriorSampler,
                                logPriorDensity = toyLogPriorDensity,
                                smcParameters$numSamples,smcParameters$numDistrbtns,
                                imageData[i,j],likeSigma  = likesigma,
                                priorSigma = priorsigma, dimParameter = 1,mu = 5)
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
          MHChainOfEstimates[[i,j,1]] = list(modelOrder  = 5,
                                             normalisingConstant = smcOuput$normalisingConstant,
                                             particles = smcOuput$particles[,smcParameters$numDistrbtns,],
                                             logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          smcOuput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                                priorSampler = toyBioPriorSampler,
                                logPriorDensity = toyLogPriorDensity,
                                smcParameters$numSamples,smcParameters$numDistrbtns,
                                imageData[i,j],likeSigma  = likesigma,
                                priorSigma = priorsigma, dimParameter = 1,mu = -5)

          currMarginalLikelihoodMatrix[i,j] =smcOuput$normalisingConstant
          MHChainOfEstimates[[i,j,1]] = list(modelOrder  = -5,
                                             normalisingConstant = smcOuput$normalisingConstant,
                                             particles = smcOuput$particles[,smcParameters$numDistrbtns,],
                                             logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }

      }
    }
    print(currMarginalLikelihoodMatrix)
    # stop("stopped")
    # print(singleMarginalEstimateMatrix[,,1])
    # print(singleMarginalEstimateMatrix[,,2])

    for(k in 2:numberOfIterations){
      # cat("\n Currently at:",k,"/",numberOfIterations)
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
          smcOutput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                                 priorSampler = toyBioPriorSampler,
                                 logPriorDensity = toyLogPriorDensity,
                                 smcParameters$numSamples,smcParameters$numDistrbtns,
                                 imageData[i,j],likeSigma  = likesigma,
                                 priorSigma = priorsigma, dimParameter = 1,mu = proposedNodeModelOrder)
          proposedNodeMarginalLikeli = smcOutput$normalisingConstant

          # cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
          # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
          # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
          # "currPottsConfigPrior: ",currPottsConfigPrior)
          # stop("stopped")
          accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
          # cat("\n acceptance Ration:",accRatio)
          if(accRatio >= runif(1)){
            # if(DEBUG2){cat("\n Accepted")}
            currModelOrderMatrix[i, j] = proposedNodeModelOrder
            currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
            MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
                                               normalisingConstant = smcOuput$normalisingConstant,
                                               particles = smcOuput$particles[,smcParameters$numDistrbtns,],
                                               logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
          }else{
            #Rejected
            MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
            # if(DEBUG2){cat("\n Rejected")}
          }
        }
      }
      MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    }
    list(MHChainOfConfigurations = MHChainOfConfigurations, MHChainOfEstimates = MHChainOfEstimates)
  }

  # n=m=3
  # J=0.5
  #
  # priorSigma = 5
  # likeSigma = 1
  #
  #
  # groundTruth =  matrix(NA,n,m)
  # groundTruth[,] = 5
  # groundTruth[2,3] = -5
  #
  # testImage = imageGeneratorToy(groundTruth)

  # testImage
#
#   nodeWisePMHSamplerToyRcppOLDTest <- function(numberOfIterations,
#                                                J = 0.6,
#                                                pottsStateSpace,
#                                                smcParameters = list(numSamples = 200,
#                                                                     numDistrbtns = 500),
#                                                priorsigma = 5, likesigma = 1, ranSeed=6
#   ){
#     DEBUG2 = F
#     DEBUG = F
#
#     dataMatrix = matrix(c(-1.7766208,  4.447216, 11.6658525
#                           ,-0.9468702, 12.697735, -8.7294234
#                           , 3.5709090,  2.305461, -0.2528843),3,3, T)
#
#     set.seed(ranSeed)
#     nodeWisePMHSamplerToyRcppOLD(dataMatrix,
#                                  numberOfIterations,
#                                  J,
#                                  pottsStateSpace,
#                                  smcParameters,
#                                  priorsigma, likesigma
#     )
#
#   }



  nodeWisePMHSamplerToyRcppTest <- function(numberOfIterations,
                                            J = 0.6,
                                            pottsStateSpace,
                                            smcParameters = list(numSamples = 200,
                                                                 numDistrbtns = 500),
                                            priorSigma = 5, likeSigma = 1, ranSeed = 6){
    DEBUG2 = F
    DEBUG = F

    dataMatrix = matrix(c(-1.7766208,  4.447216, 11.6658525
                          ,-0.9468702, 12.697735, -8.7294234
                          , 3.5709090,  2.305461, -0.2528843),3,3, T)
    # cat("dataMatrix")
    # print(dataMatrix)

    set.seed(ranSeed)
    nodeWisePMHSamplerToyRcpp(imageData = dataMatrix, numberOfIterations = numberOfIterations, J=J, pottsStateSpace = pottsStateSpace,smcParameters = smcParameters, priorSigma = priorSigma, likeSigma = likeSigma)
    # nodeWisePMHSamplerToyRcpp(imageData = dataMatrix,
    #                           numberOfIterations = numberOfIterations,
    #                           J=J,
    #                           pottsStateSpace = pottsStateSpace,
    #                           smcParameters = smcParameters,
    #                           priorsigma = priorsigma, likesigma = likesigma)$MHChainOfConfigurations
  }


  ###Actual test but takes really long, use quick one below instead for generic testin
  # expect_equal(    nodeWiseMarginalModelToy(nodeWisePMHSamplerToyRcppTest(numberOfIterations = 5, pottsStateSpace = c(-5,5),
  #                                                                     smcParameters = list(numSamples = 192,
  #                                                                                          numDistrbtns = 50),J = 0.5,
  #                                                                     priorSigma = 5, likeSigma = 1)$MHChainOfConfigurations),
  #   nodeWiseMarginalModelToy(nodeWisePMHSamplerToyRcppOLDTest(5,pottsStateSpace = c(-5,5),
  #                                                             smcParameters = list(numSamples = 192,
  #                                                                                  numDistrbtns = 50),J = 0.5,
  #                                                             priorsigma = 5, likesigma = 1)$MHChainOfConfigurations)
  # )

  ### dirty but quick test

  expect_equal(    nodeWiseMarginalModelToy(nodeWisePMHSamplerToyRcppTest(numberOfIterations = 5, pottsStateSpace = c(-5,5),
                                                                          smcParameters = list(numSamples = 192,
                                                                                               numDistrbtns = 50),J = 0.5,
                                                                          priorSigma = 5, likeSigma = 1,)$MHChainOfConfigurations),
                   matrix(c( 5  , 5 ,   5,
                             5   , 5,   -5,
                             5   , 5 ,  -5),3,3,T)
  )

})
#
#
#
# test_that("pet NWPM works", {
#
#   nodeWisePMHSamplerPETRcppOLD <- function(imageData,
#                                         numberOfIterations,
#                                         J = 0.6,
#                                         pottsStateSpace = 0,
#                                         smcParameters = list(numSamples = 200,
#                                                              numDistrbtns = 500)
#   ){
#     cat("\n smcParameters:", smcParameters$numSamples)
#     n = dim(imageData)[1]
#     m = dim(imageData)[2]
#
#     cat("n",n,"m:",m)
#
#     currModelOrderMatrix = isingibbs2Pet(niter = 10,n,m,J=J)
#     currMarginalLikelihoodMatrix = matrix(NA,n,m)
#
#     MHChainOfConfigurations = array(NA, c(n,m,numberOfIterations))
#     MHChainOfConfigurations[,,1] = currModelOrderMatrix
#
#     MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
#     # singleMarginalEstimateMatrix = array(NA,c(n,m,2))
#
#
#     cat("\n Computing marginal likelihoods using smc...")
#
#     # cat("\n Marginal Likelihoods computed....\n")
#
#     for(i in 1:n){
#       for(j in 1:m){
#
#         cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")
#         if(MHChainOfConfigurations[i,j,1] == 1){
#
#           smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,3)
#           currMarginalLikelihoodMatrix[i,j] = smcOuput$normalisingConstant
#           MHChainOfEstimates[[i,j,1]] = list(modelOrder  = 1,
#                                              normalisingConstant = smcOuput$normalisingConstant,
#                                              particles = smcOuput$particles[,smcParameters$numDistrbtns,],
#                                              logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
#         }else{
#           smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,5)
#
#           currMarginalLikelihoodMatrix[i,j] =smcOuput$normalisingConstant
#           MHChainOfEstimates[[i,j,1]] = list(modelOrder  = 2,
#                                              normalisingConstant = smcOuput$normalisingConstant,
#                                              particles = smcOuput$particles[,smcParameters$numDistrbtns,],
#                                              logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
#         }
#
#       }
#     }
#     print(currMarginalLikelihoodMatrix)
#     # print(singleMarginalEstimateMatrix[,,1])
#     # print(singleMarginalEstimateMatrix[,,2])
#
#     for(k in 2:numberOfIterations){
#       cat("\n Currently at:",k,"/",numberOfIterations)
#       for(i in 1:n){
#         for(j in 1:m){
#           #Compute the priors for the current step and the proposed step
#           # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
#           currPottsConfigPrior = pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
#           proposedNodeModelOrder = if(currModelOrderMatrix[i,j] == 1){2}else{1}
#           #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
#           # proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
#           proposedPottsConfigPrior = pottsPriorSimplified(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
#
#
#           # -------------------------------------------------------------------------
#
#           #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates
#           # currMargLik = mean(exp(unNormLogWeights[i,j,modelMatrix[i,j], ]))
#           # propMargLik  = mean(exp(unNormLogWeights[i,j,propModelOrd, ] ))
#
#           # currMargLik = dnorm(modelMatrix[i,j],0,priorSigma)*dnorm(dataMatrix[i,j],modelMatrix[i,j],likeSigma)+rnorm(1,0,1)
#           #propMargLik = dnorm(propModelOrd,0,priorSigma)*dnorm(dataMatrix[i,j],propModelOrd,likeSigma)+rnorm(1,0,1)
#
#           # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
#           smcOutput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(proposedNodeModelOrder*2)+1)
#           proposedNodeMarginalLikeli = smcOutput$normalisingConstant
#
#           cat("proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
#               "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
#               # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j],
#               "currPottsConfigPrior: ",currPottsConfigPrior)
#           accRatio = (proposedNodeMarginalLikeli * proposedPottsConfigPrior) / (currMarginalLikelihoodMatrix[i, j] * currPottsConfigPrior)
#           # cat("\n acceptance Ration:",accRatio)
#           if(accRatio >= runif(1)){
#             if(DEBUG2){cat("\n Accepted")}
#             currModelOrderMatrix[i, j] = proposedNodeModelOrder
#             currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
#             MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
#                                                normalisingConstant = smcOuput$normalisingConstant,
#                                                particles = smcOuput$particles[,smcParameters$numDistrbtns,],
#                                                logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
#           }else{
#             #Rejected
#             MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
#             if(DEBUG2){cat("\n Rejected")}
#           }
#         }
#       }
#       MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
#     }
#     list(MHChainOfConfigurations = MHChainOfConfigurations, MHChainOfEstimates = MHChainOfEstimates)
#   }
#
#
#
#
#
# })
