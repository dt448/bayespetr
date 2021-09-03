test_that("vSMC_RCpp Prior5 (+CESS adaptive x1) works", {
# Important: Test may fail due to different distribution libraries -- typically adding the flags:
# -DVSMC_HAS_CXX11LIB_FUNCTIONAL=0 -DVSMC_HAS_CXX11LIB_RANDOM=0
#  to "PKG_CXXFLAGS =" in Makevars
# should give the values quoted in the comments below.
dataAsRVec = c(-0.00501983422577359, 0.0190255021101724, -0.0661881961229502, 0.097336673064335
              ,4.68184710268953, 8.87302909558736, 9.22773434858563, 11.45993602980807, 12.3797992053381
              ,9.34326073356752, 8.81563539736783, 9.5735976815741, 8.93530528847354, 9.6095323669387
              ,9.8682761807142, 7.58235524591176, 8.36247300036969, 7.8112107791523, 8.98735650871226
              ,8.40872535117855, 6.69814564498163, 8.01129423382567, 7.37430988778437, 5.47442656474305
              ,7.06067777256831, 5.09982976254094, 5.635048519047, 4.3518222529957, 4.03829028199849
              ,3.25365377951844, 2.65614907945676, 2.002713292454)
  test_alog_config = list(mcmc_type = "MCMC",annealing_scheme = "Prior5")
  test_model_config = list(time_intervals = c(27.5, 32.5, 10.0, 10.0, 20.0, 30.0,
                                              30.0, 30.0, 30.0, 30.0, 30.0, 75.0,
                                              120.0, 120.0, 120.0, 120.0, 120.0,
                                              120.0, 120.0, 120.0, 120.0, 120.0,
                                              120.0, 210.0, 300.0, 300.0, 300.0,
                                              300.0, 300.0, 450.0, 600.0, 600.0),
                           prior_ranges = c(1e-3, 1e-5, 1e-5,
                                            0.01, 0.1, 0.1,
                                            1e-4, 1e-5, 1e-5,
                                            1e-3, 1e-1, 1e-1,
                                            1e-1,
                                            10,
                                            0,
                                            5e-1),
                           t_distribution = F,
                           decay_lowest_rate = 0)
  load_pet_conv_table("/home/den/Desktop/dataFolder/pet_conv.data")
  expect_equal(smcSamplerRcpp_vSMC(192,500,dataAsRVec,3,test_model_config,test_alog_config,1)$normalizingConstant, -57.0503574828418)
  #-56.9163521616309 - using Boost libraries and -57.47774116990816395401 on TP (-57.0503574828418 R.math)
  expect_equal(smcSamplerRcpp_vSMC(192,500,dataAsRVec,5,test_model_config,test_alog_config,231)$normalizingConstant,-44.0900556627322)
  #-50.6688638718729 - using Boost libraries and -44.84501226535390117078 on TP (-44.0900556627322 R.math)
  expect_equal(smcSamplerRcpp_vSMC(192,500,dataAsRVec,7,test_model_config,test_alog_config,546)$normalizingConstant,-49.1960220943087)
  #-59.2216015508679 - using Boost libraries adn -53.39727146860901996206 on TP (-49.1960220943087 R.math)
  # -- to check other look up table and CESS and adaptive MCMC
  load_pet_conv_table("/home/den/Desktop/dataFolder/h01031_pet_conv.data")

  test_algo_config = list(mcmc_type = "AMCMC",
                          mcmc_iters = 4,
                          annealing_scheme = "CESS",
                          cess_threshold = 0.999)
  test_model_config = list(time_intervals = c(45,  10,  10,  10,  30,  30,  30,  30,  30,  30,  30, 120, 120, 120,
                                              120, 120, 120, 120, 120, 120, 120, 120, 120, 300, 300, 300, 300, 300,
                                              300, 600, 600, 600),
                           prior_ranges = c(1e-5, 1e-5, 1e-5,
                                            1e-2, 1e-2, 1e-2,
                                            2e-4, 1e-5, 1e-5,
                                            1e-2, 1e-2, 1e-2,
                                            1e-1,
                                            10,
                                            0,
                                            5e-1),
                           t_distribution = T,
                           decay_lowest_rate = 0.0005668525)
  expect_equal(smcSamplerRcpp_vSMC(numSamples = 200,0,datum = dataAsRVec,dimParameter = 7,test_model_config,test_algo_config,  randSeed = 10)$normalizingConstant,-45.973275211228965986)
  # -96.2394859210389 - using Boost libraries and -96.0391013372735 on TP and -95.8889038276476 on MAC -96.0416729298431 on rmath #-45.97328 #Rnorm Decay fixed
})

test_that("vSMC_RCpp whole image works", {
  # Important: Test may fail due to different distribution libraries -- typically adding the flags:
  # -DVSMC_HAS_CXX11LIB_FUNCTIONAL=0 -DVSMC_HAS_CXX11LIB_RANDOM=0
  #  to "PKG_CXXFLAGS =" in Makevars
  # should give the values quoated in the comments below.

  #!!!!! IMPORTANT NEED TO LOAD IMAGE FOR THIS TO WORK!!!
  # whole image experiments can be found in 210405_EXP1 directory
  # load_pet_conv_table("/home/den/Desktop/dataFolder/h01031_pet_conv.data")0

  # test_algo_config = list(mcmc_type = "AMCMC",
  #                         mcmc_iters = 4,
  #                         annealing_scheme = "CESS")
  # test_model_config = list(time_intervals = c(45,  10,  10,  10,  30,  30,  30,  30,  30,  30,  30, 120, 120, 120,
  #                                             120, 120, 120, 120, 120, 120, 120, 120, 120, 300, 300, 300, 300, 300,
  #                                             300, 600, 600, 600),
  #                          prior_ranges = c(1e-5, 1e-5, 1e-5,
  #                                           1e-2, 1e-2, 1e-2,
  #                                           7e-4, 1e-5, 1e-5,
  #                                           1e-2, 1e-2, 1e-2,
  #                                           1e-1,
  #                                           10,
  #                                           0,
  #                                           5e-1),
  #                          t_distribution = T,
  #                          decay_lowest_rate = 0.0005668525)
  # set.seed(6)
  # expect_equal(smcSamplerRcpp_vSMC(numSamples = 1000, cess_threshold = 0.999,datum = h01032_4950[1,1,],dimParameter = 3,test_model_config,test_algo_config),77.7546614576012)
  #77.7546614576012 - using Boost libraries and 77.7409783232295 on TP and 77.7195261415226 on MAC
})



