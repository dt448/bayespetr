test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("smc_toy normalising constant works", {
  smcTest <- function(seedFixed = 6){

    priorSigma = 100
    likeSigma = 10
    mu=0

    dataSet = 5.00614

    set.seed(seedFixed)
    # smcSamplerToy(numSamples = 192,numDistrbtns = 50, dataSet, likeSigma, priorSigma, varDimsn = 1,mu)$normalizingConstant
    smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
                          logPriorDensity = toyLogPriorDensity,
                          priorSampler = toyBioPriorSampler,
                          numSamples = 192,
                          numDistributions = 50,
                          dimParameter = 1,
                          datum = dataSet,
                          priorSigma = priorSigma,
                          likeSigma = likeSigma,
                          mu = 0
  )$normalisingConstant

  }
  expect_equal(smcTest(), 0.004135367502)
  expect_equal(smcTest(60), 0.003863509152)
  expect_equal(smcTest(1543), 0.004015785737)
})




test_that("smc_toy_rcpp normalising constant works", {
  smcTest <- function(seedFixed = 6){

    priorSigma = 100
    likeSigma = 10
    mu=0

    dataSet = 4.99503

    set.seed(seedFixed)
    test = smcSamplerToyRcpp(192,50,dataSet,1,priorSigma,likeSigma,mu)$normalisingConstant
  }
  expect_equal(smcTest(), 0.004136001)
  expect_equal(smcTest(60), 0.003855626)
  # expect_equal(smcTest(125), 0.00403637)
  expect_equal(smcTest(1543), 0.004034848)
})



test_that("smc_toy sample works", {
  smcTest <- function(seedFixed = 6, samplePosition = 60, distributionPosition = 16){

    priorSigma = 100
    likeSigma = 10
    mu=0

    dataSet = 5.00614

    set.seed(seedFixed)
    # smcSamplerToy(numSamples = 192,numDistrbtns = 50, dataSet, likeSigma, priorSigma, varDimsn = 1,mu)$samples[samplePosition,distributionPosition,1]
    smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
               logPriorDensity = toyLogPriorDensity,
               priorSampler = toyBioPriorSampler,
               numSamples = 192,
               numDistributions = 50,
               dimParameter = 1,
               datum = dataSet,
               priorSigma = priorSigma,
               likeSigma = likeSigma,
               mu = 0
    )$particles[samplePosition,distributionPosition,1]
  }
  expect_equal(smcTest(), 70.48114143)
  expect_equal(smcTest(seedFixed = 1453), 3.411542705)
  expect_equal(smcTest(seedFixed = 60,samplePosition = 100), -212.4660762)
  expect_equal(smcTest(seedFixed = 423, distributionPosition = 23), -20.76194613)
})

test_that("smc_toy weights works", {
  smcTest <- function(seedFixed = 61, samplePosition = 60, distributionPosition = 23, weightType = "barLogWeights"){

    priorSigma = 100
    likeSigma = 10
    mu=0

    dataSet = 5.00614

    set.seed(seedFixed)
    # smcOutput = smcSamplerToy(numSamples = 192,numDistrbtns = 50, dataSet, likeSigma, priorSigma, varDimsn = 1,mu)
    smcOutput = smcSampler(logLikelihoodDensity = toyLogLikelihoodDensity,
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
    if (weightType == "barLogWeights"){
      smcOutput$barLogWeights[samplePosition,distributionPosition]
    }else if (weightType == "incremWeights"){
     smcOutput$incremWeights[samplePosition,distributionPosition]
    }

  }
  expect_equal(smcTest(), -5.333944539)
  expect_equal(smcTest( weightType = "incremWeights"), -0.01877593566)
  expect_equal(smcTest(seedFixed = 213, samplePosition = 32), -5.095686227)
})

test_that("smc_pet_Rcpp normalising constant works", {

  smcTest <- function(seedFixed = 6, dimParameter = 5,numDistributions = 500){
    phi = c(0.001752992, 0.004947008)
    theta = c(0.0111990653, 0.0005009347)

    dataSet = c (0.0008199750,0.0017871899 , 0.0018674304, -0.0006880768, -0.0004680085,  3.4993148899,  7.0149475981,  8.3965827114 , 9.0373437567 ,
                 9.2884931257 ,9.3662051640  ,9.3268443765  ,9.4911005701,  9.1717982210,  8.8833974168,  8.5500954157,  8.3578397516,  8.1274822345,
                 7.6957669549 , 7.4726012098  ,7.1896656682 , 6.8339554226,  6.5436026896  ,6.0752194184,  5.2872657085,  4.7246544526  ,4.1290334444,
                 3.6299445843  ,3.1518110475,  2.6224344293,  1.9760028477,  1.5049115747)

    timeFrames = c(0.0,27.5,60.0,70.0,80.0,100.0,130.0,160.0,190.0,220.0,250.0,280.0,355.0,
               475.0,595.0,715.0,835.0,955.0,1075.0,1195.0,1315.0,1435.0,1555.0,1675.0,1885.0,2185.0,
               2485.0,2785.0,3085.0,3385.0,3835.0,4435.0,5035.0)

    set.seed(seedFixed)
    # log(smcSampler(200,500,y = dataSet, timeFrames = timeFrames, modelId)$normalisingConstant)
    log(smcSamplerRcpp(200,numDistributions, dataSet , timeFrames, dimParameter = dimParameter)$normalisingConstant)
  }
  #33.02146 when using RCpp
  expect_equal(smcTest(), 39.0039157)
  expect_equal(smcTest(dimParameter = 7,numDistributions = 50), -13.43542661)
  expect_equal(smcTest(dimParameter = 3,numDistributions = 100), -3.96957648)
})

test_that("smc_pet_R normalising constant works", {
  ## Rcpp implementation makes this obselete but can be run if needed,
  ## WARNING: Will take long time to run.
  smcTest <- function(seedFixed = 6, dimParameter = 5,numDistributions = 500){
    phi = c(0.001752992, 0.004947008)
    theta = c(0.0111990653, 0.0005009347)

    dataSet = c (0.0008199750,0.0017871899 , 0.0018674304, -0.0006880768, -0.0004680085,  3.4993148899,  7.0149475981,  8.3965827114 , 9.0373437567 ,
                 9.2884931257 ,9.3662051640  ,9.3268443765  ,9.4911005701,  9.1717982210,  8.8833974168,  8.5500954157,  8.3578397516,  8.1274822345,
                 7.6957669549 , 7.4726012098  ,7.1896656682 , 6.8339554226,  6.5436026896  ,6.0752194184,  5.2872657085,  4.7246544526  ,4.1290334444,
                 3.6299445843  ,3.1518110475,  2.6224344293,  1.9760028477,  1.5049115747)

    set.seed(seedFixed)
    # log(smcSampler(200,500,y = dataSet, timeFrames = timeFrames, modelId)$normalisingConstant)

    log(smcSampler(petLogLikelihoodDensity,petLogBioPriorDensity,petBioPriorSampler,200,numDistributions,dimParameter = dimParameter,datum = dataSet,timeFrames = timeFrames)$normalisingConstant)
  }
  #33.02146 when using RCpp
  # expect_equal(smcTest(), 39.0039157)
  # expect_equal(smcTest(dimParameter = 7,numDistributions = 50), -13.43542661)
  # expect_equal(smcTest(dimParameter = 3,numDistributions = 100), -3.96957648)

})
