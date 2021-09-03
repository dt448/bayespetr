######## Toy Model Demo  for R/Rcpp smc sampler--###############################
################################################################################

### Preamble

## Run this before starting code -- used for debugging to removed before release
DEBUG=F
DEBUG2 = T


############################### Simulated Data Generation ######################
################################################################################

# First we fix the parameters and the observed data

# Set the hyperparameters:

priorSigma = 100
likeSigma = 10
mu =0

### Data generated from generating distribution
#dataSet = rnorm(1,mu = 0,likeSigma)

# Fix observed data
dataSet = 5.00614

############################ R SMC Sampler #####################################
################################################################################

# Fix random seeed
set.seed(6)

# Run sampler
test2 = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                          logPriorDensity = toyLogPriorDensity,
                          priorSampler = toyBioPriorSampler,
                          numSamples = 192,
                          numDistributions = 50,
                          dimParameter = 1,
                          datum = dataSet,
                          priorSigma = priorSigma,
                          likeSigma = likeSigma,
                          mu = 0
)

# Call normalising constant
test2$normalisingConstant
########## Everything is running correctly if the above returns :0.004135368



############################ RCppp #############################################
################################################################################

#### Next to check the Rcpp version

# First Fix random seeed
set.seed(6)

# Then run sampler
rcppTest = smcSamplerToyRcpp(numSamples = 192,
                             numDistributions = 50,
                             datum = dataSet,
                             dimParameter = 1,
                             priorSigma = priorSigma,
                             likeSigma = likeSigma,
                             mu = mu)

# Call normalising constant
rcppTest$normalisingConstant
########## As before you should get the same number :0.004135368

######################### Computing The True normalising Constant###############
################################################################################

# This function will return the true normalising constant, since it can be compu
# ted in this toy example

bayespetr::toyModelTrueNormalizingConst(dataSet,
                                        modelOrder = 0,
                                        priorSigma = priorSigma,
                                        likelihoodSigma = likeSigma)



######-- Check Asymptotics--###################################################

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


######--Checking True Value--##################################################

toyModelTrueNormalizingConst(dataSet,modelOrder = 0,priorSigma = priorSigma,likelihoodSigma = likeSigma)

sigmaFG = (priorSigma^2*likeSigma^2)/(priorSigma^2+likeSigma^2)
sigmaFG
scaleFactor = exp(0.5*(
  ((priorSigma^2*dataSet^2)/(likeSigma^2*(likeSigma^2+priorSigma^2))) - ((dataSet^2)/(likeSigma^2))

))

scaleFactor

trueValue = (sqrt(2*pi*sigmaFG)*scaleFactor)/(sqrt(2*pi*likeSigma^2)*sqrt(2*pi*priorSigma^2))
trueValue
