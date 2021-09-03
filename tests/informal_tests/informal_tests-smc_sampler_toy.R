# Preamble
DEBUG=F
DEBUG2 = T

priorSigma = 100
likeSigma = 10

dataSet = 5.00614

## To check convergence as we increase the number of particles.
estimate2 = rep(NA,50)
for(k in 1:50){
  test2 = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                            logPriorDensity = toyLogPriorDensity,
                            priorSampler = toyBioPriorSampler,
                            numSamples = (50+k^2),
                            numDistributions = 100,
                            dimParameter = 1,
                            datum = dataSet,
                            priorSigma = priorSigma,
                            likeSigma = likeSigma,
                            mu = 0
  )

  estimate2[k] = test2$normalisingConstant
  print(k)
  # print(test)
}

realValue = bayespetr::toyModelTrueNormalizingConst(dataSet,
                                        modelOrder = 0,
                                        priorSigma = priorSigma,
                                        likelihoodSigma = likeSigma)

plot(estimate2[1:50],
     xlab = "Number of particles used (50+k^2)",
     main = "Toy model normalising constant: Increasing number of particles,
     True value shown in red")
abline(h=realValue, col="red")
mean(estimate2)
abline(h=mean(estimate2),lty=2)


# To check the cumulative mean with more estimates.

estimateMean = rep(NA,200)
for(k in 1:200){
  test2 = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                     logPriorDensity = toyLogPriorDensity,
                     priorSampler = toyBioPriorSampler,
                     numSamples = 200,
                     numDistributions = 50,
                     dimParameter = 1,
                     datum = dataSet,
                     priorSigma = priorSigma,
                     likeSigma = likeSigma,
                     mu = 0
  )

  estimateMean[k] = test2$normalisingConstant
  print(k)
  # print(test)
}

plot(cumsum(estimateMean)/seq(1,200), xlab = "Number of estimations",
     main = " Toy model normalising constant: Cumulative averages of normalising
     constant estimates, True Value in red")
abline(h=realValue, col = "red")
mean(estimateMean)
log(var(estimateMean))

plot(density(estimateMean),
     main = "Denisty of estimates of normalising constnat, realvalue = 0.003964,
     mean = 0.003968,
     log variance = -17.03569855")
abline(v=realValue, col = "red")
