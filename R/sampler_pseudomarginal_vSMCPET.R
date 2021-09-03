load_table_env <- new.env(parent = emptyenv())

load_table_env$VMSC_PETCONV_TABLE_LOADED = F

#'@export
load_pet_conv_table <- function(pet_conv_table_dir){
  if(!file.exists(pet_conv_table_dir)){
    stop("File not found!")
  }
  cpp_load_conv_table(pet_conv_table_dir);
  old <- load_table_env$VMSC_PETCONV_TABLE_LOADED
  load_table_env$VMSC_PETCONV_TABLE_LOADED <- T
  invisible(old)
}

#'@export
smcSamplerRcpp_vSMC <-function(numSamples = 200,
                               prior5_iter_num = 500,
                               datum,
                               dimParameter,
                               config_model = list(time_intervals,
                                                   prior_ranges = c(1e-3, 1e-5, 1e-5,
                                                                    0.01, 0.1, 0.1,
                                                                    1e-4, 1e-5, 1e-5,
                                                                    1e-3, 1e-1, 1e-1,
                                                                    1e-1,
                                                                    10,
                                                                    0,
                                                                    5e-1),
                                                   t_distribution = T,
                                                   decay_lowest_rate),
                               config_algo = list( mcmc_type = "AMCMC",
                                                   mcmc_iters=1,
                                                   annealing_scheme = "Prior5",
                                                   cess_threshold),
                               randSeed = NULL) {
  # Sanitizing inputs
  if(!load_table_env$VMSC_PETCONV_TABLE_LOADED == T){
    stop("Please load pet conv look up table first.")
  }
  if(is.null(config_algo$mcmc_type) | is.null(config_algo$annealing_scheme) | is.null(config_model$time_intervals)){
    stop("Error: Invalid config")
  }
  if(is.null(randSeed)){
    randSeed =  floor(runif(1, min=1, max=1000))
  }else{
    set.seed(randSeed)
  }
  if(is.null(config_algo$mcmc_iters)){
    config_algo$mcmc_iters = 1
  }
  if(!config_algo$annealing_scheme == "User"&!config_algo$annealing_scheme == "CESS"&!config_algo$annealing_scheme == "Prior5"){
    stop("Pease give correct annealing scheme: 'User', 'CESS' or 'Prior5'")
  }
  if(config_algo$annealing_scheme == "User"){
    if(is.null(config_algo$user_annealing_scheme))
      stop("Pease input an annealing scheme vector")
    if(config_algo$user_annealing_scheme[length(config_algo$user_annealing_scheme)] < 1)
      stop("Last alpha value should be 1 or greater")
  }
  modelOrder = (dimParameter-1)/2
  config_algo$t_distribution = config_model$t_distribution

  #passing data and other model information to vSMC C++
  dataFromR(datum,config_model)

  vSMCRCppSampler(numSamples, prior5_iter_num, modelOrder, R_config = config_algo, randSeed);
  # smcSamplerCpp(numSamples, numDistributions, datum, timeFrames, dimParameter)
}


#'@export
pet_smc_indep <- function(imageData,
                          smcParameters = list(numSamples = 200,
                                               numDistrbtns = 500),
                          config_model = list(time_intervals,
                                              prior_ranges = c(1e-3, 1e-5, 1e-5,
                                                               0.01, 0.1, 0.1,
                                                               1e-4, 1e-5, 1e-5,
                                                               1e-3, 1e-1, 1e-1,
                                                               1e-1,
                                                               10,
                                                               0,
                                                               5e-1),
                                              t_distribution = T,
                                              decay_lowest_rate),
                          config_algo = list( mcmc_type = "AMCMC",
                                              mcmc_iters=1,
                                              annealing_scheme = "Prior5",
                                              cess_threshold)
                          ){

  cat("\n smcParameters:", smcParameters$numSamples)
  n = dim(imageData)[1]
  m = dim(imageData)[2]
  selectedModel = matrix(NA,n,m)
  V_D_output = matrix(NA,n,m)
  cat("n:", n, "m:", m)

  for(i in 1:n){
    cat("\n Currently at, i:",i,"/",n)
    for(j in 1:m){
      # smcOutput_1 = smcSamplerToyRcpp(numSamples = smcParameters$numSamples,
      #                                 numDistributions = smcParameters$numDistrbtns,
      #                                 datum = imageData[i,j],
      #                                 dimParameter = 1,
      #                                 priorSigma = priorSigma,
      #                                 likeSigma = likeSigma,
      #                                 mu = 5)
      smcOutput_1 = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                     smcParameters$numDistrbtns,
                                     imageData[i,j,],
                                     3,
                                     config_model,
                                     config_algo)
      smcOutput_2 = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                        smcParameters$numDistrbtns,
                                        imageData[i,j,],
                                        5,
                                        config_model,
                                        config_algo)
      smcOutput_3 = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                        smcParameters$numDistrbtns,
                                        imageData[i,j,],
                                        7,
                                        config_model,
                                        config_algo)
      if(smcOutput_1$normalizingConstant>smcOutput_2$normalizingConstant & smcOutput_1$normalizingConstant>smcOutput_3$normalizingConstant){
        selectedModel[i,j] = 1
        V_D_output[i,j] = smcOutput_1$volumeOfDistribution
      }else if (smcOutput_2$normalizingConstant>smcOutput_3$normalizingConstant){
        selectedModel[i,j] = 2
        V_D_output[i,j] = smcOutput_2$volumeOfDistribution
      }else{
        selectedModel[i,j] = 3
        V_D_output[i,j] = smcOutput_3$volumeOfDistribution
      }
    }
  }

  list(selectedModel = selectedModel, V_D = V_D_output)
}




#'@export
NWPM_SamplerPETRcpp_vSMC <- function(imageData,
                                      numberOfIterations,
                                      J = 0.6,
                                      pottsStateSpace = c(1,2,3),
                                      smcParameters = list(numSamples = 200,
                                                           numDistrbtns = 500),
                                     config_model = list(time_intervals,
                                                         prior_ranges = c(1e-3, 1e-5, 1e-5,
                                                                          0.01, 0.1, 0.1,
                                                                          1e-4, 1e-5, 1e-5,
                                                                          1e-3, 1e-1, 1e-1,
                                                                          1e-1,
                                                                          10,
                                                                          0,
                                                                          5e-1),
                                                         t_distribution = T,
                                                         decay_lowest_rate),
                                     config_algo = list( mcmc_type = "AMCMC",
                                                         mcmc_iters=1,
                                                         annealing_scheme = "Prior5",
                                                         cess_threshold),
                                     initial_state_info = NULL,
                                     SAVEDIR = NULL){

  # Sanitizing inputs
  if(!load_table_env$VMSC_PETCONV_TABLE_LOADED == T){
    stop("Please load pet conv look up table first.")
  }

  if(is.null(config_algo$mcmc_type) | is.null(config_algo$annealing_scheme) | is.null(config_model$time_intervals)){
    stop("Error: Invalid config")
  }

  n = dim(imageData)[1]
  m = dim(imageData)[2]

  currModelOrderMatrix = matrix(NA,n,m)
  # this matrix will keep record of current likelihood -- used for calculating
  # acceptance ratio
  currMarginalLikelihoodMatrix = matrix(NA,n,m)
  numStates = length(pottsStateSpace)

  MHChainOfConfigurations = V_D_Chain = array(NA, c(n,m,numberOfIterations))
  # MHChainOfConfigurations[,,1] = currModelOrderMatrix

  # MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
  # singleMarginalEstimateMatrix = array(NA,c(n,m,2))


  # cat("\n Computing marginal likelihoods using smc...")

  # cat("\n Marginal Likelihoods computed....\n")
  cat("\n Using Potts state space:", pottsStateSpace)
  #if initalisation states,likelihoods etc is NOT provided, sample from Potts model
  if(is.null(initial_state_info)){
    cat("\n Initializing states... \n")
    pb = txtProgressBar(1, n, style=3)
    #sample MODEL ORDERS from Potts using Gibbs Sampler
    currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,
                      n = n,
                      m = m,
                      J=J,
                      stateSpace = pottsStateSpace)
    #Record the initial state.
    MHChainOfConfigurations[,,1] = currModelOrderMatrix
    # Next compute the marginal likeihoods of these model orders
    for(i in 1:n){
      for(j in 1:m){

        # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

        # smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
        #Compute marginal and V_D using currModelOrder (inital state)
        smcOuput = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                       smcParameters$numDistrbtns,
                                       imageData[i,j,],
                                       (MHChainOfConfigurations[i,j,1]*2)+1,
                                       config_model,
                                       config_algo)
        # save to currMarg
        currMarginalLikelihoodMatrix[i,j] = smcOuput$normalizingConstant
        V_D_Chain[i,j,1] = smcOuput$volumeOfDistribution
        # MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
        #                                    normalisingConstant = smcOuput$normalisingConstant,
        #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
        #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
      }
      setTxtProgressBar(pb, i)
    }
    cat("\n Initializing complete... \n")
    # print(currMarginalLikelihoodMatrix)
  }else{
    # If intial states, marginals are provided.
    #sanatizing inputs
    if(is.null(initial_state_info$model_orders)
       | is.null(initial_state_info$marginal_estimates)
       | is.null(initial_state_info$marginal_estimates)
    ){
      stop("Invalid initial state information.")
    }
    cat("\n Using user initial states")
    currModelOrderMatrix = initial_state_info$model_orders
    if(max(currModelOrderMatrix)>max(pottsStateSpace)){
      stop("Potts state space is smaller than given initial state sapce.")
    }
    #Record the initial state, set currMarg to given states.
    MHChainOfConfigurations[,,1] = currModelOrderMatrix
    currMarginalLikelihoodMatrix = initial_state_info$marginal_estimates
    V_D_Chain[,,1] = initial_state_info$volume_distributions
  }
  # print(singleMarginalEstimateMatrix[,,1])
  # print(singleMarginalEstimateMatrix[,,2])

  cat("\r Pseudo-Marginal MCMC Chain Progress ... \n")


  for(k in 2:numberOfIterations){
    cat("\r \n Currently at:",k,"/",numberOfIterations, " iterations \n")
    pb = txtProgressBar(1, n*m, style=3)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        ## We want to propose new states (i.e. other than the current state), so remove current state
        proposalStateSpace = pottsStateSpace[!pottsStateSpace == currModelOrderMatrix[i,j]]
        proposedNodeModelOrder = NA
        ## Propose from the propossalStateSpace
        if(numStates == 2){
          proposedNodeModelOrder = proposalStateSpace[1]
        }else{
          proposedNodeModelOrder = sample(proposalStateSpace,1)
        }
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        # proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # Compute prior density of proposed model order
        proposedPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        # smcOutput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(proposedNodeModelOrder*2)+1)

        # Compute the marginal likelihood and v_D of proposed model order
        smcOuput = smcSamplerRcpp_vSMC(numSamples = smcParameters$numSamples,
                            prior5_iter_num = smcParameters$numDistrbtns,
                            datum = imageData[i,j,],
                            dimParameter = (proposedNodeModelOrder*2)+1,
                            config_model = config_model,
                            config_algo = config_algo)

        proposedNodeMarginalLikeli = smcOuput$normalizingConstant
        proposedV_D = smcOuput$volumeOfDistribution

        # cat("\n proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currPottsConfigPrior: ",currPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j])

        # STOPP INCORRECT -- the marginalLikelihood is given by log in vSMC so double check this !!!!
        # accRatio = (exp(proposedNodeMarginalLikeli) * proposedPottsConfigPrior) / (exp(currMarginalLikelihoodMatrix[i, j]) * currPottsConfigPrior)
        # Compute the acceptance ratio
        accRatio = (proposedNodeMarginalLikeli+log(proposedPottsConfigPrior)) - (currMarginalLikelihoodMatrix[i, j]+log(currPottsConfigPrior))

        if(exp(accRatio) >= runif(1)){
          if(DEBUG2){cat("\n Accepted \n")}
          # if accepted update the current model order and marginal likelihood matrix
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          V_D_Chain[i,j,k] = proposedV_D
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected - keep the same order and marginal likelihood matrix
          V_D_Chain[i,j,k] = V_D_Chain[i,j,(k-1)]
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
        setTxtProgressBar(pb, ((i-1)*m)+j)
      }
    }
    #record whole image for that iteration
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    #save incremntal progress if requested
    if(!is.null(SAVEDIR)){
      saveRDS(MHChainOfConfigurations[,,1:k],paste(SAVEDIR,"MHChainOfConfigurations_iter_",k,".RDS",sep = ""))
      saveRDS(V_D_Chain[,,1:k],paste(SAVEDIR,"V_D_Chain_iter_",k,".RDS",sep = ""))
    }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations, V_D_Chain = V_D_Chain)#, MHChainOfEstimates = MHChainOfEstimates)
}

## NWPM-SE approximation for PET
#'@export
pet_NWPM_SE <- function(imageData,
                        numberOfIterations,
                        J,
                        pottsStateSpace,
                        smcParameters = list(numSamples = 200,
                                             numDistrbtns = 500),
                        config_model = list(time_intervals,
                                            prior_ranges = c(1e-3, 1e-5, 1e-5,
                                                             0.01, 0.1, 0.1,
                                                             1e-4, 1e-5, 1e-5,
                                                             1e-3, 1e-1, 1e-1,
                                                             1e-1,
                                                             10,
                                                             0,
                                                             5e-1),
                                            t_distribution = T,
                                            decay_lowest_rate),
                        config_algo = list( mcmc_type = "AMCMC",
                                            mcmc_iters=1,
                                            annealing_scheme = "Prior5",
                                            cess_threshold),
                        SAVE_DIR = NULL
){

  # Sanitizing inputs
  if(!load_table_env$VMSC_PETCONV_TABLE_LOADED == T){
    stop("Please load pet conv look up table first.")
  }

  if(is.null(config_algo$mcmc_type) | is.null(config_algo$annealing_scheme) | is.null(config_model$time_intervals)){
    stop("Error: Invalid config")
  }

  n = dim(imageData)[1]
  m = dim(imageData)[2]

  currModelOrderMatrix = matrix(NA,n,m)
  # this matrix will keep record of current likelihood -- used for calculating
  # acceptance ratio
  currMarginalLikelihoodMatrix = matrix(NA,n,m)
  numStates = length(pottsStateSpace)

  MHChainOfConfigurations = V_D_Chain = array(NA, c(n,m,numberOfIterations))
  # MHChainOfConfigurations[,,1] = currModelOrderMatrix

  # MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
  SE_marginal_estimates  = array(NA,c(n,m,numStates))
  V_D_estimates  = array(NA,c(n,m,numStates))


  # cat("\n Computing marginal likelihoods using smc...")

  # cat("\n Marginal Likelihoods computed....\n")
  cat("\n Using Potts state space:", pottsStateSpace)
  #if initalisation states,likelihoods etc is NOT provided, sample from Potts model
  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,
                                           n = n,
                                           m = m,
                                           J=J,
                                           stateSpace = pottsStateSpace)
  #Record the initial state.
  MHChainOfConfigurations[,,1] = currModelOrderMatrix
  # Next compute the marginal likeihoods of these model orders
  for(state in 1:numStates){
    for(i in 1:n){
    pb = txtProgressBar(1, n, style=3)
      for(j in 1:m){

        # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

        # smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
        #Compute marginal and V_D using currModelOrder (inital state)
        smcOuput = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                       smcParameters$numDistrbtns,
                                       imageData[i,j,],
                                       (pottsStateSpace[state]*2)+1,
                                       config_model,
                                       config_algo)
        SE_marginal_estimates[i,j,state] = smcOuput$normalizingConstant
        V_D_estimates[i,j,state] = smcOuput$volumeOfDistribution
        # save to currMarg
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalizingConstant
        }
        # MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
        #                                    normalisingConstant = smcOuput$normalisingConstant,
        #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
        #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
      }
      setTxtProgressBar(pb, i)
    }
  }
  cat("\r Pseudo-Marginal MCMC Chain Progress ... \n")


  for(k in 2:numberOfIterations){
    cat("\r \n Currently at:",k,"/",numberOfIterations, " iterations \n")
    # pb = txtProgressBar(1, n*m, style=3)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        ## We want to propose new states (i.e. other than the current state), so remove current state
        proposalStateSpace = pottsStateSpace[!pottsStateSpace == currModelOrderMatrix[i,j]]
        proposedNodeModelOrder = NA
        ## Propose from the propossalStateSpace
        if(numStates == 2){
          proposedNodeModelOrder = proposalStateSpace[1]
        }else{
          proposedNodeModelOrder = sample(proposalStateSpace,1)
        }
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        # proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # Compute prior density of proposed model order
        proposedPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        # smcOutput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(proposedNodeModelOrder*2)+1)

        # Compute the marginal likelihood and v_D of proposed model order
        # smcOuput = smcSamplerRcpp_vSMC(numSamples = smcParameters$numSamples,
        #                                prior5_iter_num = smcParameters$numDistrbtns,
        #                                datum = imageData[i,j,],
        #                                dimParameter = (proposedNodeModelOrder*2)+1,
        #                                config_model = config_model,
        #                                config_algo = config_algo)

        proposedNodeMarginalLikeli = SE_marginal_estimates[i,j,which(proposedNodeModelOrder==pottsStateSpace)]
        # proposedV_D = smcOuput$volumeOfDistribution

        # cat("\n proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currPottsConfigPrior: ",currPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j])

        # STOPP INCORRECT -- the marginalLikelihood is given by log in vSMC so double check this !!!!
        # accRatio = (exp(proposedNodeMarginalLikeli) * proposedPottsConfigPrior) / (exp(currMarginalLikelihoodMatrix[i, j]) * currPottsConfigPrior)
        # Compute the acceptance ratio
        accRatio = (proposedNodeMarginalLikeli+log(proposedPottsConfigPrior)) - (currMarginalLikelihoodMatrix[i, j]+log(currPottsConfigPrior))

        if(exp(accRatio) >= runif(1)){
          if(DEBUG2){cat("\n Accepted \n")}
          # if accepted update the current model order and marginal likelihood matrix
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          # V_D_Chain[i,j,k] = proposedV_D
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected - keep the same order and marginal likelihood matrix
          # V_D_Chain[i,j,k] = V_D_Chain[i,j,(k-1)]
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
        # setTxtProgressBar(pb, ((i-1)*m)+j)
      }
    }
    #record whole image for that iteration
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    if(!is.null(SAVE_DIR) && k%%1000 == 0){
      saveRDS(MHChainOfConfigurations,paste(SAVE_DIR,"MHChainOfConfigurations_PROGRESS.RDS",sep = ""))
      # saveRDS(V_D_Chain,paste(SAVE_DIR,"V_D_Chain_PROGRESS.RDS",sep = ""))
    }
    #save incremntal progress if requested
  #   if(!is.null(SAVEDIR)){
  #     saveRDS(MHChainOfConfigurations[,,1:k],paste(SAVEDIR,"MHChainOfConfigurations_iter_",k,".RDS",sep = ""))
  #     saveRDS(V_D_Chain[,,1:k],paste(SAVEDIR,"V_D_Chain_iter_",k,".RDS",sep = ""))
  #   }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations,SE_marginal_estimates = SE_marginal_estimates, V_D_estimates = V_D_estimates)#, MHChainOfEstimates = MHChainOfEstimates)

}

## NWPM-MA approximation for PET
#'@export
pet_NWPM_MA <- function(imageData,
                        numberOfIterations,
                        J,
                        pottsStateSpace,
                        smcParameters = list(numSamples = 200,
                                             numDistrbtns = 500),
                        config_model = list(time_intervals,
                                            prior_ranges = c(1e-3, 1e-5, 1e-5,
                                                             0.01, 0.1, 0.1,
                                                             1e-4, 1e-5, 1e-5,
                                                             1e-3, 1e-1, 1e-1,
                                                             1e-1,
                                                             10,
                                                             0,
                                                             5e-1),
                                            t_distribution = T,
                                            decay_lowest_rate),
                        config_algo = list( mcmc_type = "AMCMC",
                                            mcmc_iters=1,
                                            annealing_scheme = "Prior5",
                                            cess_threshold),
                        REFRESH=F,refreshConst
){

  # Sanitizing inputs
  if(!load_table_env$VMSC_PETCONV_TABLE_LOADED == T){
    stop("Please load pet conv look up table first.")
  }

  if(is.null(config_algo$mcmc_type) | is.null(config_algo$annealing_scheme) | is.null(config_model$time_intervals)){
    stop("Error: Invalid config")
  }

  n = dim(imageData)[1]
  m = dim(imageData)[2]

  currModelOrderMatrix = matrix(NA,n,m)
  # this matrix will keep record of current likelihood -- used for calculating
  # acceptance ratio
  currMarginalLikelihoodMatrix = matrix(NA,n,m)
  numStates = length(pottsStateSpace)

  MHChainOfConfigurations = V_D_Chain = array(NA, c(n,m,numberOfIterations))
  # MHChainOfConfigurations[,,1] = currModelOrderMatrix

  # MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
  SE_marginal_estimates  = array(NA,c(n,m,numStates))
  V_D_estimates  = array(NA,c(n,m,numStates))


  # cat("\n Computing marginal likelihoods using smc...")

  # cat("\n Marginal Likelihoods computed....\n")
  cat("\n Using Potts state space:", pottsStateSpace)
  #if initalisation states,likelihoods etc is NOT provided, sample from Potts model
  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,
                                           n = n,
                                           m = m,
                                           J=J,
                                           stateSpace = pottsStateSpace)
  #Record the initial state.
  MHChainOfConfigurations[,,1] = currModelOrderMatrix
  # Next compute the marginal likeihoods of these model orders
  for(state in 1:numStates){
    for(i in 1:n){
      pb = txtProgressBar(1, n, style=3)
      for(j in 1:m){

        # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

        # smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
        #Compute marginal and V_D using currModelOrder (inital state)
        smcOuput = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                       smcParameters$numDistrbtns,
                                       imageData[i,j,],
                                       (pottsStateSpace[state]*2)+1,
                                       config_model,
                                       config_algo)
        SE_marginal_estimates[i,j,state] = smcOuput$normalizingConstant
        V_D_estimates[i,j,state] = smcOuput$volumeOfDistribution
        # save to currMarg
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalizingConstant
        }
        # MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
        #                                    normalisingConstant = smcOuput$normalisingConstant,
        #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
        #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
      }
      setTxtProgressBar(pb, i)
    }
  }
  cat("\r Pseudo-Marginal MCMC Chain Progress ... \n")
  print(currMarginalLikelihoodMatrix)
  # print(SE_marginal_estimates)


  for(k in 2:numberOfIterations){
    cat("\r \n Currently at:",k,"/",numberOfIterations, " iterations \n")
    pb = txtProgressBar(1, n*m, style=3)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        ## We want to propose new states (i.e. other than the current state), so remove current state
        proposalStateSpace = pottsStateSpace[!pottsStateSpace == currModelOrderMatrix[i,j]]
        proposedNodeModelOrder = NA
        ## Propose from the propossalStateSpace
        if(numStates == 2){
          proposedNodeModelOrder = proposalStateSpace[1]
        }else{
          proposedNodeModelOrder = sample(proposalStateSpace,1)
        }
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        # proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # Compute prior density of proposed model order
        proposedPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        # smcOutput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(proposedNodeModelOrder*2)+1)

        # Compute the marginal likelihood and v_D of proposed model order
        # smcOuput = smcSamplerRcpp_vSMC(numSamples = smcParameters$numSamples,
        #                                prior5_iter_num = smcParameters$numDistrbtns,
        #                                datum = imageData[i,j,],
        #                                dimParameter = (proposedNodeModelOrder*2)+1,
        #                                config_model = config_model,
        #                                config_algo = config_algo)

        proposedNodeMarginalLikeli = SE_marginal_estimates[i,j,which(proposedNodeModelOrder==pottsStateSpace)]
        # proposedV_D = smcOuput$volumeOfDistribution

        # cat("\n proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currPottsConfigPrior: ",currPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j])

        # STOPP INCORRECT -- the marginalLikelihood is given by log in vSMC so double check this !!!!
        # accRatio = (exp(proposedNodeMarginalLikeli) * proposedPottsConfigPrior) / (exp(currMarginalLikelihoodMatrix[i, j]) * currPottsConfigPrior)
        # Compute the acceptance ratio
        accRatio = (proposedNodeMarginalLikeli+log(proposedPottsConfigPrior)) - (currMarginalLikelihoodMatrix[i, j]+log(currPottsConfigPrior))

        if(exp(accRatio) >= runif(1)){
          if(DEBUG2){cat("\n Accepted \n")}
          # if accepted update the current model order and marginal likelihood matrix
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          V_D_Chain[i,j,k] = V_D_estimates[i,j,proposedNodeModelOrder]
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected - keep the same order and marginal likelihood matrix
          V_D_Chain[i,j,k] = V_D_Chain[i,j,(k-1)]
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
        setTxtProgressBar(pb, ((i-1)*m)+j)
      }
    }
    #record whole image for that iteration
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    #save incremntal progress if requested
    #   if(!is.null(SAVEDIR)){
    #     saveRDS(MHChainOfConfigurations[,,1:k],paste(SAVEDIR,"MHChainOfConfigurations_iter_",k,".RDS",sep = ""))
    #     saveRDS(V_D_Chain[,,1:k],paste(SAVEDIR,"V_D_Chain_iter_",k,".RDS",sep = ""))
    #   }
    if(REFRESH == T && k %% refreshConst == 0 && k!=numberOfIterations){
    for(state in 1:numStates){
      for(i in 1:n){
        pb = txtProgressBar(1, n, style=3)
        for(j in 1:m){

          # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

          # smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
          #Compute marginal and V_D using currModelOrder (inital state)
          smcOuput = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                         smcParameters$numDistrbtns,
                                         imageData[i,j,],
                                         (pottsStateSpace[state]*2)+1,
                                         config_model,
                                         config_algo)
          SE_marginal_estimates[i,j,state] = smcOuput$normalizingConstant
          V_D_estimates[i,j,state] = smcOuput$volumeOfDistribution
          # save to currMarg
          if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
            currMarginalLikelihoodMatrix[i,j] = smcOuput$normalizingConstant
          }
          # MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }
        setTxtProgressBar(pb, i)
      }
    }
    }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations,SE_marginal_estimates = SE_marginal_estimates, V_D_estimates = V_D_estimates)#, MHChainOfEstimates = MHChainOfEstimates)

}


## NWPM-MA approximation for PET
#'@export
pet_NWPM_MA2 <- function(imageData,
                        numberOfIterations,
                        J,
                        pottsStateSpace,
                        smcParameters = list(numSamples = 200,
                                             numDistrbtns = 500),
                        config_model = list(time_intervals,
                                            prior_ranges = c(1e-3, 1e-5, 1e-5,
                                                             0.01, 0.1, 0.1,
                                                             1e-4, 1e-5, 1e-5,
                                                             1e-3, 1e-1, 1e-1,
                                                             1e-1,
                                                             10,
                                                             0,
                                                             5e-1),
                                            t_distribution = T,
                                            decay_lowest_rate),
                        config_algo = list( mcmc_type = "AMCMC",
                                            mcmc_iters=1,
                                            annealing_scheme = "Prior5",
                                            cess_threshold),
                        REFRESH=F,refreshConst,
                        SAVE_DIR = NULL
){

  # Sanitizing inputs
  if(!load_table_env$VMSC_PETCONV_TABLE_LOADED == T){
    stop("Please load pet conv look up table first.")
  }

  if(is.null(config_algo$mcmc_type) | is.null(config_algo$annealing_scheme) | is.null(config_model$time_intervals)){
    stop("Error: Invalid config")
  }

  n = dim(imageData)[1]
  m = dim(imageData)[2]

  currModelOrderMatrix = matrix(NA,n,m)
  # this matrix will keep record of current likelihood -- used for calculating
  # acceptance ratio
  currMarginalLikelihoodMatrix = matrix(NA,n,m)
  numStates = length(pottsStateSpace)

  MHChainOfConfigurations = V_D_Chain = array(NA, c(n,m,numberOfIterations))
  # MHChainOfConfigurations[,,1] = currModelOrderMatrix

  # MHChainOfEstimates = array(list(NULL), c(n,m,numberOfIterations))
  SE_marginal_estimates  = array(NA,c(n,m,numStates))
  V_D_estimates  = array(NA,c(n,m,numStates))


  # cat("\n Computing marginal likelihoods using smc...")

  # cat("\n Marginal Likelihoods computed....\n")
  cat("\n Using Potts state space:", pottsStateSpace)
  #if initalisation states,likelihoods etc is NOT provided, sample from Potts model
  currModelOrderMatrix = pottsGibbsSampler(numIteration = 10,
                                           n = n,
                                           m = m,
                                           J=J,
                                           stateSpace = pottsStateSpace)
  #Record the initial state.
  MHChainOfConfigurations[,,1] = currModelOrderMatrix
  # Next compute the marginal likeihoods of these model orders
  for(state in 1:numStates){
    for(i in 1:n){
      pb = txtProgressBar(1, n, style=3)
      for(j in 1:m){

        # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

        # smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
        #Compute marginal and V_D using currModelOrder (inital state)
        smcOuput = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                       smcParameters$numDistrbtns,
                                       imageData[i,j,],
                                       (pottsStateSpace[state]*2)+1,
                                       config_model,
                                       config_algo)
        SE_marginal_estimates[i,j,state] = smcOuput$normalizingConstant
        V_D_estimates[i,j,state] = smcOuput$volumeOfDistribution
        # save to currMarg
        if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
          currMarginalLikelihoodMatrix[i,j] = smcOuput$normalizingConstant
        }
        # MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
        #                                    normalisingConstant = smcOuput$normalisingConstant,
        #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
        #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
      }
      setTxtProgressBar(pb, i)
    }
  }
  cat("\r Pseudo-Marginal MCMC Chain Progress ... \n")
  print(currMarginalLikelihoodMatrix)
  # print(SE_marginal_estimates)


  for(k in 2:numberOfIterations){
    cat("\r \n Currently at:",k,"/",numberOfIterations, " iterations \n")
    # pb = txtProgressBar(1, n*m, style=3)
    for(i in 1:n){
      for(j in 1:m){
        #Compute the priors for the current step and the proposed step
        # currPottsConfigPrior = isingPrior(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)
        currPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,currModelOrderMatrix[i,j],currModelOrderMatrix,J)

        ## We want to propose new states (i.e. other than the current state), so remove current state
        proposalStateSpace = pottsStateSpace[!pottsStateSpace == currModelOrderMatrix[i,j]]
        proposedNodeModelOrder = NA
        ## Propose from the propossalStateSpace
        if(numStates == 2){
          proposedNodeModelOrder = proposalStateSpace[1]
        }else{
          proposedNodeModelOrder = sample(proposalStateSpace,1)
        }
        #propModelOrd = sample(modelOrders[modelOrders!=modelMatrix[i,j]],1,prob =c(0.5,0.5)) #Generate random proposed step.
        # proposedPottsConfigPrior = isingPrior(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # Compute prior density of proposed model order
        proposedPottsConfigPrior = bayespetr:::pottsPriorSimplified(i,j,proposedNodeModelOrder,currModelOrderMatrix,J)
        # -------------------------------------------------------------------------

        #Look up and Assign the marginal likelihood uisng the (stored) SMC estimates

        # currMargLik  =  smcCESS(numSamples,numDistrbtns,dataSet = dataMatrix[i,j],likeSigma,priorSigma = ,varDimsn = 1,mu = modelMatrix[i,j])
        # smcOutput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(proposedNodeModelOrder*2)+1)

        # Compute the marginal likelihood and v_D of proposed model order
        # smcOuput = smcSamplerRcpp_vSMC(numSamples = smcParameters$numSamples,
        #                                prior5_iter_num = smcParameters$numDistrbtns,
        #                                datum = imageData[i,j,],
        #                                dimParameter = (proposedNodeModelOrder*2)+1,
        #                                config_model = config_model,
        #                                config_algo = config_algo)

        proposedNodeMarginalLikeli = SE_marginal_estimates[i,j,which(proposedNodeModelOrder==pottsStateSpace)]
        # proposedV_D = smcOuput$volumeOfDistribution

        # cat("\n proposedNodeMarginalLikelii: ",proposedNodeMarginalLikeli,
        # "proposedPottsConfigPrior: ",proposedPottsConfigPrior,
        # "currPottsConfigPrior: ",currPottsConfigPrior,
        # "currMarginalLikelihoodMatrix: ",currMarginalLikelihoodMatrix[i,j])

        # STOPP INCORRECT -- the marginalLikelihood is given by log in vSMC so double check this !!!!
        # accRatio = (exp(proposedNodeMarginalLikeli) * proposedPottsConfigPrior) / (exp(currMarginalLikelihoodMatrix[i, j]) * currPottsConfigPrior)
        # Compute the acceptance ratio
        accRatio = (proposedNodeMarginalLikeli+log(proposedPottsConfigPrior)) - (currMarginalLikelihoodMatrix[i, j]+log(currPottsConfigPrior))

        if(exp(accRatio) >= runif(1)){
          if(DEBUG2){cat("\n Accepted \n")}
          # if accepted update the current model order and marginal likelihood matrix
          currModelOrderMatrix[i, j] = proposedNodeModelOrder
          currMarginalLikelihoodMatrix[i, j] = proposedNodeMarginalLikeli
          V_D_Chain[i,j,k] = V_D_estimates[i,j,proposedNodeModelOrder]
          # MHChainOfEstimates[[i,j,k]] = list(modelOrder  = proposedNodeModelOrder,
          #                                    normalisingConstant = smcOuput$normalisingConstant,
          #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
          #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
        }else{
          #Rejected - keep the same order and marginal likelihood matrix
          V_D_Chain[i,j,k] = V_D_Chain[i,j,(k-1)]
          # MHChainOfEstimates[[i,j,k]] = MHChainOfEstimates[[i,j,k-1]]
          if(DEBUG2){cat("\n Rejected")}
        }
        # setTxtProgressBar(pb, ((i-1)*m)+j)
      }
    }
    #record whole image for that iteration
    MHChainOfConfigurations[ , ,k] = currModelOrderMatrix
    #save incremntal progress if requested
    if(REFRESH == T && k %% refreshConst == 0 && k!=numberOfIterations){
      if(!is.null(SAVE_DIR)){
        saveRDS(MHChainOfConfigurations,paste(SAVE_DIR,"MHChainOfConfigurations_PROGRESS.RDS",sep = ""))
        saveRDS(V_D_Chain,paste(SAVE_DIR,"V_D_Chain_PROGRESS.RDS",sep = ""))
      }
    for(state in 1:numStates){
      for(i in 1:n){
        pb = txtProgressBar(1, n, style=3)
        for(j in 1:m){

          # cat("\n Currently at:",j+(m*(i-1)),"/",n*m, "\n")

          # smcOuput = smcSamplerRcpp(smcParameters$numSamples,smcParameters$numDistrbtns,imageData[i,j,],timeFrames,(MHChainOfConfigurations[i,j,1]*2)+1)
          #Compute marginal and V_D using currModelOrder (inital state)
          smcOuput = smcSamplerRcpp_vSMC(smcParameters$numSamples,
                                         smcParameters$numDistrbtns,
                                         imageData[i,j,],
                                         (pottsStateSpace[state]*2)+1,
                                         config_model,
                                         config_algo)
          accratio = exp(smcOuput$normalizingConstant-SE_marginal_estimates[i,j,state])
         if(accratio>= runif(1)){
           SE_marginal_estimates[i,j,state] = smcOuput$normalizingConstant
           V_D_estimates[i,j,state] = smcOuput$volumeOfDistribution
           # save to currMarg
           if(MHChainOfConfigurations[i,j,1] == pottsStateSpace[state]){
             currMarginalLikelihoodMatrix[i,j] = smcOuput$normalizingConstant
           }
           # MHChainOfEstimates[[i,j,1]] = list(modelOrder  = MHChainOfConfigurations[i,j,1],
           #                                    normalisingConstant = smcOuput$normalisingConstant,
           #                                    particles = smcOuput$particles[,smcParameters$numDistrbtns,],
           #                                    logWeights = smcOuput$normLogWeights[,smcParameters$numDistrbtns])
         }
        }
        setTxtProgressBar(pb, i)
      }
    }
    }
  }
  list(MHChainOfConfigurations = MHChainOfConfigurations,SE_marginal_estimates = SE_marginal_estimates, V_D_Chain = V_D_Chain,V_D_estimates = V_D_estimates)#, MHChainOfEstimates = MHChainOfEstimates)

}



#--- V_D

vd_posterior <- function(vd_nwpm_chain,nwpm_chain,slected_mo, selec_dim = dim(slected_mo)){
  wid = selec_dim[1]
  len = selec_dim[2]
  output = matrix(NA,wid,len)
  for(i in 1:wid){
    for(j in 1:len){

      output[i,j] = mean(vd_nwpm_chain[i,j,nwpm_chain[i,j,] == slected_mo[i,j]])

    }
  }
  output
}

vd_mod_average <-function(vd_nwpm_chain,selec_dim = dim(vd_nwpm_chain)){
  wid = selec_dim[1]
  len = selec_dim[2]
  output = matrix(NA,wid,len)

  for(i in 1:wid){
    output[i,] =  rowMeans(vd_nwpm_chain[i,,])
  }
  output
}
