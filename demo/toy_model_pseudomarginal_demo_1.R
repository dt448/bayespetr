n=m=20
J=0.5

priorSigma = 5
likeSigma = 1


groundTruth =  matrix(NA,n,m)
groundTruth[,] = 5
groundTruth[5:10,1:4] = -5

groundTruth

groundTruth =  matrix(NA,n,m)
groundTruth[,] = 5
groundTruth[15:19,2:5] = -5

groundTruth[3,17] = -5
groundTruth[4,16:18]  = -5
groundTruth[5,15:19]  = -5
groundTruth[6,14:20]  = -5
groundTruth[7,15:19]  = -5
groundTruth[8,16:18]  = -5
groundTruth[9,17] = -5

groundTruth[18:20,16:19] = -5
groundTruth[17,18:19] = -5
groundTruth[16,16:19] = -5

groundTruth

testImage = imageGeneratorToy(groundTruth)


# testFullMarginal = nodeWisePMHSampler(testImage$dataMatrix, numberOfIterations = 100, J=0.5, smcParameters = list(numSamples = 200,  numDistrbtns=15),toyModelParameters = list(likeSigma = likeSigma, priorSigma = priorSigma))


chain1 = nodeWisePMHSamplerToyRcpp(imageData = testImage$dataMatrix,
                                   numberOfIterations = 50, J=5.0,
                                   pottsStateSpace = c(-5,5),
                                   smcParameters = list(numSamples = 20,
                                                        numDistrbtns = 40),
                                   priorSigma = priorSigma, likeSigma = likeSigma)


DEBUG2 = F
DEBUG = F
nodeWiseMarginalModelToy(chain1$MHChainOfConfigurations)

correctPercentage(nodeWiseMarginalModelToy(chain1$MHChainOfConfigurations),groundTruth)


test2 = nodeWisePMHApproxSampler(testImage$dataMatrix, numberOfIterations = 500, J=0.7, smcParameters = list(numSamples = 200,  numDistrbtns=15),toyModelParameters = list(likeSigma = likeSigma, priorSigma = priorSigma))

DEBUG2 = T
DEBUG = F
nodeWiseMarginalModelToy(test2$MHChainOfConfigurations)
nodeWiseMarginalModelToy(testFullMarginal)

test2[,,3]


correctPercentage(nodeWiseMarginalModelToy(test2),groundTruth)
