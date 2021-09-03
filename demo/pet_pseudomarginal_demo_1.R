# Preamble
DEBUG=F
DEBUG2 = T


timeFrames = c(0.0,27.5, 60.0, 70.0, 80.0, 100.0, 130.0, 160.0, 190.0, 220.0, 250.0,280.0,
               355.0, 475.0,595.0,715.0,835.0,955.0,1075.0,1195.0,1315.0,1435.0,
               1555.0,1675.0,1885.0,2185.0, 2485.0,2785.0,3085.0,3385.0,3835.0,
               4435.0, 5035.0)


n=m=20

#=============================== Ground truth 1 ================================
# Use this for smaller ground truth image: for faster computation time
# used mostly for informal testing....

n=m=5
groundTruth =  matrix(NA,n,m)
groundTruth[,] = 2
groundTruth[3:4,3:4] = 1

#=============================== Ground truth 2 ================================
groundTruth =  matrix(NA,n,m)
groundTruth[,] = 1
groundTruth[15:19,2:5] = 2

groundTruth[3,17] = 2
groundTruth[4,16:18]  = 2
groundTruth[5,15:19]  = 2
groundTruth[6,14:20]  = 2
groundTruth[7,15:19]  = 2
groundTruth[8,16:18]  = 2
groundTruth[9,17] = 2

groundTruth[18:20,16:19] = 2
groundTruth[17,18:19] = 2
groundTruth[16,16:19] = 2

#============================= pseud-marginal algorithm=========================

groundTruth


petImageData = imageGeneratorPet(groundTruth)

saveRDS(petImageData,"/home/den/Documents/Simulation Experiments/040420EXP1 PETPottsExperimentFullvsSE/petImageData20by20.RDS")

petImageData[20,20,]

chain = nodeWisePMHApproxSamplerPET(petImageData,
                                    200,
                                    smcParameters = list(numSamples = 200,
                                                         numDistrbtns = 500),J = 0.5)

chain2 = nodeWisePMHSamplerPET(petImageData,
                                    50,
                                    smcParameters = list(numSamples = 200,
                                                         numDistrbtns = 500),J = 0.5)

chain4 = nodeWisePMHSamplerPETRcpp(petImageData,
                               50,
                               smcParameters = list(numSamples = 140,
                                                    numDistrbtns = 400),J = 0.5)
SEEstimates = nodeWisePMHApproxSamplerPETRcpp(petImageData,100,J=0.6)

load_pet_conv_table("/home/den/Desktop/dataFolder/h01031_pet_conv.data")

demo_algo_config = list(mcmc_type = "AMCMC",
                       mcmc_iters = 4,
                       annealing_scheme = "Prior5")

demo_model_config = list(time_intervals = c(45,  10,  10,  10,  30,  30,  30,  30,  30,  30,  30, 120, 120, 120,
                                           120, 120, 120, 120, 120, 120, 120, 120, 120, 300, 300, 300, 300, 300,
                                           300, 600, 600, 600),
                        prior_ranges = c(1e-5, 1e-5, 1e-5,
                                         1e-2, 1e-2, 1e-2,
                                         7e-4, 1e-5, 1e-5,
                                         1e-2, 1e-2, 1e-2,
                                         1e-1,
                                         10,
                                         0,
                                         5e-1),
                        t_distribution = T,
                        decay_lowest_rate = 0.0005668525)

vSMCChain = NWPM_SamplerPETRcpp_vSMC(imageData = petImageData,
                                     numberOfIterations = 5,
                                     smcParameters = list(numSamples = 140,
                                                          numDistrbtns = 400),
                                     J = 0.5,
                                     config_model= demo_model_config,
                                     config_algo = demo_algo_config)

# initial_state_info$model_orders
# currMarginalLikelihoodMatrix = initial_state_info$marginal_estimates
# V_D_Chain[,,1] = initial_state_info$volume_distributions

demo_initial_state = list(model_orders = matrix(1,n,m),
                          marginal_estimates = matrix(120,n,m),
                          volume_distributions = matrix(10,n,m))

vSMCChain = NWPM_SamplerPETRcpp_vSMC(imageData = petImageData,
                                     numberOfIterations = 5,
                                     smcParameters = list(numSamples = 140,
                                                          numDistrbtns = 400),
                                     J = 0.5,
                                     config_model= demo_model_config,
                                     config_algo = demo_algo_config,
                                     initial_state_info = demo_initial_state)


DEBUG = F
DEBUG2 = F
OPTIMCONV = T

nodeWiseMarginalModelPet(chain$MHChainOfConfigurations)
nodeWiseMarginalModelPet(chain2$MHChainOfConfigurations)
nodeWiseMarginalModelPet(chain3$MHChainOfConfigurations)
nodeWiseMarginalModelPet(chain4$MHChainOfConfigurations)
nodeWiseMarginalModelPet(SEEstimates$MHChainOfConfigurations)
