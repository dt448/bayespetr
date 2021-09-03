##-----------------PREAMPLE--------------------------------------------------
OPTIMCONV = T
DEBUG = F
DEBUG2 = T

# Alternative parameters to use if needed-----
# phi=c(0.0044406807, 0.0001010392 )
# theta = c(0.0004518293, 0.0027728739)

# This is what is presented in
# phi = c(0.001752992, 0.004947008)
# theta = c(0.0111990653, 0.0005009347)

#This reordering as required by the biologically informed priors.
phi = c(0.004947008, 0.001752992)
theta = c(0.0005009347,0.0111990653)

timeFrames = c(0.0,27.5,60.0,70.0,80.0,100.0,130.0,160.0,190.0,220.0,250.0,280.0,355.0,
               475.0,595.0,715.0,835.0,955.0,1075.0,1195.0,1315.0,1435.0,1555.0,1675.0,1885.0,2185.0,
               2485.0,2785.0,3085.0,3385.0,3835.0,4435.0,5035.0)

#Simulated data
data1 = c (0.0008199750,0.0017871899 , 0.0018674304, -0.0006880768, -0.0004680085,  3.4993148899,  7.0149475981,  8.3965827114 , 9.0373437567 ,
          9.2884931257 ,9.3662051640  ,9.3268443765  ,9.4911005701,  9.1717982210,  8.8833974168,  8.5500954157,  8.3578397516,  8.1274822345,
          7.6957669549 , 7.4726012098  ,7.1896656682 , 6.8339554226,  6.5436026896  ,6.0752194184,  5.2872657085,  4.7246544526  ,4.1290334444,
          3.6299445843  ,3.1518110475,  2.6224344293,  1.9760028477,  1.5049115747)


######--Demo of R SMC Sampler--##################################################
# Here we will demonstrate how to use the (R) smc sampler for PET
# Begin by fixing the seed

set.seed(6)

#  Next call the smc sampler function (note the parameters),
# in this example we will use  2-compartmental model for PET data with
# with biologically informed priors.
# The 2-compartment model results in a 5-dimensional( 2*num of comp + 1)
# parameter space. Hence, the dimParater=5, the rest of the arguements should
# be self-explanatory -- see ?smcSampler for more information.

test = smcSampler(
  logLikelihoodDensity = petLogLikelihoodDensity,
  logPriorDensity = petLogBioPriorDensity,
  priorSampler = petBioPriorSampler,
  numSamples = 200,
  numDistributions = 500,
  dimParameter = 5,
  datum = data1,
  timeFrames = timeFrames
)

# If everything runs correctly you should get that:
# log(test$normalisingConstant) = 39.00392

log(test$normalisingConstant)

######--Demo for RCpp SMC Sampler--############################################
# We will now repeat the above for the RCpp implementation --- we should get the
# the same normalizing constant

set.seed(6)
testRcpp = smcSamplerRcpp(200,500,data1,timeFrames,5)
log(testRcpp$normalisingConstant)

## To call the (R)Cpp function directly we run the following
set.seed(6)
# We will use R to generate the intial set of particles, and then pass these
# onto the Cpp sampler --- the above call (line 56) using the function
# smcSamplerRcpp does this automatically.
rTakeit(petBioPriorSampler(5,200,timeFrames))
testCpp = smcSamplerCpp(200,500,data1,timeFrames,5)
log(testCpp$normalisingConstant)

# If everything runs correctly you should get that:
# log(testRccp$normalisingConstant) = 39.00392


# To check Zhou et al paper, run the following
# set.seed(6)
# phi=c(0.0044406807, 0.0001010392 )
# theta = c(0.0004518293, 0.0027728739)
#
# data1Yan = dataSimltr(timeFrames,phi,theta,0.64)$Y
#
# rTakeit(petBioPriorSampler(5,200,timeFrames))
# testRccp = smcSamplerRcpp(200,500,dataSet = data1Yan,timeFrames,5)
# log(testRccp$normalisingConstant)
## log(testRccp$normalisingConstant) = -33.83422
