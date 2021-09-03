test

y = dataSimltr(timeFrames,phi,theta,0.64)

test = smcSamplerGeneral(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,500,dimParameter = 5,datum = y$Y,timeFrames = timeFrames)
test1 = smcSamplerGeneral(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,500,dimParameter = 3,datum = y$Y,timeFrames = timeFrames)

test$normalisingConstant>test1$normalisingConstant

testR = smcSamplerGeneral(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,500,dimParameter = 5,datum = y$Y,timeFrames = timeFrames)
testR1 = smcSamplerGeneral(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,500,dimParameter = 3,datum = y$Y,timeFrames = timeFrames)

testR$normalisingConstant>testR1$normalisingConstant

log(testR$normalisingConstant)
log(testR1$normalisingConstant)
plot(testR1$ESS)
plot(colMeans(testR1$acceptanceRates[,2:500]))

rTakeit(petBioPriorSampler(5,300,timeFrames))

    testRc2 = smcSamplerRcpp(200,500,y$Y,timeFrames,5)
    testRc1 = smcSamplerRcpp(200,500,y$Y,timeFrames,3)

    testRc2$normalisingConstant>testRc1$normalisingConstant

numTrials = 30
estimates1Comp = matrix(NA,numTrials,4)
estimates2Comp = matrix(NA,numTrials,4)
estimates3Comp = matrix(NA,numTrials,4)
noiseVector = c(0.01,0.64,1.28,5.12)
noiseVector = c(0.01,0.04,0.08,0.16)
# noiseVector = c(0.64)

for(noiseLeves in noiseVector){

  y = dataSimltr(timeFrames,phi,theta,noiseLeves)


  for(i in 1:numTrials){
    # testRc3 = smcSamplerRcpp(200,500,y$Y,timeFrames,7)
    testRc2 = smcSamplerRcpp(200,500,y$Y,timeFrames,5)
    testRc1 = smcSamplerRcpp(200,500,y$Y,timeFrames,3)
    # estimates3Comp[i,which(noiseLeves==noiseVector) ] = testRc3$normalisingConstant
    estimates2Comp[i,which(noiseLeves==noiseVector) ] = testRc2$normalisingConstant
    estimates1Comp[i,which(noiseLeves==noiseVector) ] = testRc1$normalisingConstant
    saveRDS(estimates1Comp,"~/Desktop/Estimate1newSeed.RDS")
    saveRDS(estimates2Comp,"~/Desktop/Estimate2newSeed.RDS")
    # saveRDS(estimates3Comp,"~/Desktop/Estimate3Compe.RDS")
  }
}

for(noiseLeves in noiseVector){

  y = dataSimltr(timeFrames,phi,theta,noiseLeves)


  for(i in 1:numTrials){
    # testRc3 = smcSamplerRcpp(200,500,y$Y,timeFrames,7)
    testRc2 = smcSamplerGeneral(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,50,5,y$Y,timeFrames=timeFrames)
    testRc1 = smcSamplerGeneral(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,500,3,y$Y,timeFrames=timeFrames)
    # estimates3Comp[i,which(noiseLeves==noiseVector) ] = testRc3$normalisingConstant
    estimates2Comp[i,which(noiseLeves==noiseVector) ] = testRc2$normalisingConstant
    estimates1Comp[i,which(noiseLeves==noiseVector) ] = testRc1$normalisingConstant
  }
  saveRDS(estimates1Comp,"~/Desktop/Estimate1CompeRSMC.RDS")
  saveRDS(estimates2Comp,"~/Desktop/Estimate2CompeRSMC.RDS")
  saveRDS(estimates3Comp,"~/Desktop/Estimate3CompeRSMC.RDS")
}
var(log(estimates1Comp))
var(log(estimates2Comp))
mean(log(estimates2Comp))
mean(log(estimates1Comp))

plot(testRc$ESS)
plot(testRc1$ESS)
plot(testRc1$acceptanceRates)
plot(testRc$acceptanceRates)

testRc$particles[,500,5]
testRc1$particles[,500,3]

log(testRc$normalisingConstant)
log(testRc1$normalisingConstant)

estimatesSmc = colSums(test$samples[,500,]*exp(test$normLogWeights[,500]))
sum(estimatesSmc[1:2]/estimatesSmc[3:4])

sum(phi/theta)


noiseLevels = c(0.01,0.16,0.64,1.28,5.12)


for( n in noiseLevels){

  simulatedData = dataSimltr(timeFrames,phi,theta,noiseLevel = n)
  # simulatedData = readRDS("~/Desktop/variabilityOfSingleEstimateTest/simData1.RDS")
  noiseFreeData = simulatedData$C_T
  ySMC = simulatedData$Y

  saveRDS(simulatedData,file =paste("~/Desktop/SimulationExperiments/PETMargLikeliHoodVarianceExperiment/1000Dis960Particles/simData"
                                    ,toString(n),".RDS",sep = ""))

}


estimates1000Dis960Part = matrix(NA,30,5)
colnames(estimates500Dis192Part) = c("0.01","0.16","0.64","1.28","5.12")

# toString(noiseLevels[5])

for(k in 11:20){

  for(n in noiseLevels){
    simulatedData = readRDS(paste("~/Google Drive/Project Von Nuemann/SimulationExperiments/PETMargLikeliHoodVarianceExperiment/5000Dis192Particles/simData"
                                  ,toString(n),".RDS",sep = ""))
    noiseFreeData = simulatedData$C_T
    ySMC = simulatedData$Y
    test = smcCESS(numSamples =192,numDistrbtns=5000,y = ySMC,timeFrames = timeFrames,modelId=2)

    prod1=mean(exp(test$incremWeights[,1]))
    prod1
    for(i in 2:(which(test$alphaVec == 1))){
      prod1 = prod1 * sum(exp(test$incremWeights[,i]+test$barLogWeights[,i-1]))

    }
    # estimates1000Dis960Part[k] = prod1
    saveRDS(test,file = paste("~/Google Drive/Project Von Nuemann/SimulationExperiments/PETMargLikeliHoodVarianceExperiment/5000Dis192Particles/",toString(n),"/trail",toString(k),".RDS",sep=""))
    # saveRDS(estimates1000Dis960Part,file = paste("~/Google Drive/Project Von Nuemann/SimulationExperiments/PETMargLikeliHoodVarianceExperiment/1000Dis960Particles/estimates1000Dis960Part.RDS",sep=""))

  }
}
