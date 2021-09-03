input_fun_env <- new.env(parent = emptyenv())
load("data/h01031_parentplasma_input.rda",input_fun_env)
InputFun = splinefun(input_fun_env$h01031_parentplasma_input$time, input_fun_env$h01031_parentplasma_input$CP)

# InputFun = splinefun(SimulatedPlasmaInput$time, SimulatedPlasmaInput$CP)
#plot(0:5500,InputFun(0:5500),type = 'l', xlab = "Time t (secs)", ylab = "Concentration (kBq/ml)")

BasisFun <- function(x, uptime, theta){
  exp(-theta * (uptime - x)) * InputFun(x)
}

#This is the tissue time-activity function.
C_T = function(t,phiVec,thetaVec){
  sum<-0
  for(i in 1:length(phiVec)){
    sum<- sum + abs(phiVec[i]*integrate(BasisFun, lower=0, upper=t,theta=thetaVec[i]
                                        ,uptime=t,
                                        subdivisions=1e8, stop.on.error=FALSE)$value)
    # if(DEBUG){cat('\n','C_T Summation', i ,' Element:',sum,'\n')} # prints out each element of the summation, used for debugging
  }

  # if(DEBUG){cat('\n','Final C_T Output:',sum,'\n')}
  return(sum)
}


# DEBUG = F
# ct = rep(NA,5000+1)
# ct = rep(NA,(length(TimeFrame)))
# for(i in 1:(length(TimeFrame))){
#   ct[i] = C_T(TimeFrame[i],c(0.0044406807,0.0014582801,0.0001010392),c(0.0004518293,0.0107752968, 0.0027728739))
# }
#
# plot(TimeFrame,ct,type = 'l', xlab = "Time t (secs)", ylab = "Concentration (kBq/ml)")

#Use this function to simulate noisy (and noise-free) data
dataSimltr = function(timeFrames,phi,theta,noiseLevel){

  y = vector('numeric',length(timeFrames)-1)
  ct = vector('numeric',length(timeFrames)-1)

  for(i in 1:length(y)){
    ct[i] = C_T(timeFrames[i+1],phi,theta)
  }
  deltaStar = timeFrames[(which(ct==max(ct))+1)]-timeFrames[which(ct==max(ct))]
  kStar =deltaStar/ max(ct)
  lambda = 1/(noiseLevel*kStar)

  for(i in 1:length(y)){
    y[i]= ct[i]+(sqrt(ct[i]/(timeFrames[i+1]-timeFrames[i]))*rnorm(1,0,sqrt(1/lambda)))
  }
  return (list (Y=y , C_T=ct,Lambda = lambda, B = log(lambda), kStar = kStar,TrueValues = c(phi,theta)))
}


Convolution <- function(theta, t){
  abs(integrate(BasisFun, lower=0, upper=t, uptime=t, theta=theta,
                subdivisions=1e8, stop.on.error=FALSE)$value)
}

convMatElem = function(theta,i){
  #Optimized but extremely slow due to using sysdata, for testing only.
  convVector = as.vector(t( convMat_Matrix ))
  convVector[(theta*length(timeFrames))+i]
}


ComputeConv <- function(theta_a, theta_b, theta_step, timeFrames)
{
  Theta <- seq(theta_a, theta_b, theta_step)
  m <- length(Theta)
  convMul = m
  n <- length(timeFrames)

  ConvList <- list()

  for (i in 1:n) {
    cat("Time frame: ", i, "out of ", n,"\n")
    tf <- timeFrames[i]
    ConvList[[i]] <- parallel::mclapply(Theta, Convolution, t = tf)
  }

  ConvMat <- matrix(rep(0, m * n), ncol = n)
  for (j in 1:n) {
    for (i in 1:m) {
      ConvMat[i,j] <- ConvList[[j]][[i]]
    }
  }

  list(ConvMat=ConvMat, Theta=Theta, TimeFrames = timeFrames)
}

## THIS IS OUTDATED NEED TO BE FIXED
# Lookup table generated using this code, timeFrames used are identical to simulation experiments.
# ConvMat <- ComputeConv(1e-4, 1, 1e-5,timeFrames)
# save(ConvMat,file = "~/Desktop/gitpull/mcmcimplem/ConvMat.Rdata")

# intpolFitValueRold = function(i,phi,theta){


# convMul = dim(bayespetr:::convMat_Matrix)[1]

#   #We will estimate the fitted values using linear interpolation using the convoultion grid matrix
#   fv = 0
#   modelID = length(theta)
#   for (d in 1:modelID){
#     thetaCur = theta[d] * convMul
#     thetaLow = floor(thetaCur)+1
#     # if(DEBUG){ cat("\n thetaLow:",thetaLow)}
#     thetaUp = thetaLow+1
#     # if(DEBUG){ cat("\n thetaUp:",thetaUp)}
#     lowMulti = phi[d] * (thetaCur - thetaLow)
#     uppMulti = phi[d] * (thetaUp - thetaCur)
#     # if(DEBUG){ cat("\n uppMulti",uppMulti)}
#     # if(DEBUG){ cat("\n lowMulti",lowMulti)}
#     # if(DEBUG){ cat("\n convMatElem(thetaUp,i):",convMatElem(thetaUp,i))}
#     fv = fv + (lowMulti*convMatElem(thetaUp,i) + uppMulti*convMatElem(thetaLow,i))
#     # if(DEBUG){ cat("\n fv:",fv)}
#   }
#   return(fv)
#
# }

#Bio-Informed Priors

phivPrior = function(phi)
{
  prior = 0
  if( (phi>=(10^-5)) && (phi<= 1) ){
    prior = 1
  }
  return(prior)

}

thetavPrior = function(theta)
{
  prior = 0
  if((theta>=10^-4) && (theta<= 1)){
    prior = 1
  }
  return(prior)
}

dBioPrior=function(phiVec,thetaVec){

  dim = length(phiVec)
  # for(i in 1:dim){
  #   if(phiVec[i] <0 || thetaVec[i]<0){
  #     cat("\n Phi:",phiVec)
  #     cat("\n Theta:",thetaVec)
  #     # stop("negative Proposal!")
  #   }
  # }

  #print(dim)
  phiPriors = thetaPriors = rep(0,dim)

  phiMean = 3e-3

  for(i in 1:dim){
    #cat("LoopStart")
    phiPriors[i] = truncnorm::dtruncnorm(phiVec[i],a = 1e-5,b = 1e-2,mean = phiMean,sqrt(1e-3))
    #cat("phiPriors",phiPriors[i])
    if(i==1){
      thetaPriors[i] = truncnorm::dtruncnorm(thetaVec[i],a = 2e-4,b = 1e-2,mean =(phiVec[i]/15) ,sqrt(1e-2))
      phiMean = 1e-3
      #cat("ThetaPriors",thetaPriors[i])
    }
    else if(i==2){
      thetaPriors[i] = truncnorm::dtruncnorm(thetaVec[i],a = thetaVec[1],b = 6e-2,mean =(phiVec[i]/4) ,sqrt(1e-2))

      #cat("thetaPriors",thetaPriors[i])
    }
    else if(i==3){
      thetaPriors[i] = truncnorm::dtruncnorm(thetaVec[i],a = thetaVec[i-1],b = 6e-2,mean =phiVec[i],sqrt(1e-2))
    }

  }

  prod(c(phiPriors,thetaPriors))
}

alpha<-beta<-10^-3

logLikelihood = function(data,timeFrames,phiVec,thetaVec,b){

  # We will need a summand as we are looking at the log transformed liklihood.
  logLHood = 0

  for(i in 1:length(data)){
    # We compute the tissue concentration(C_T - at the end point of each time frame)


    if(OPTIMCONV == T){
      c_t = intpolFitValue(i,phiVec,thetaVec)
      #if(c_t==0){return<--Inf}else{return<- return - (lambda/2)*((tVec[i+1]-tVec[i])/c_t)*((dataVec[i]-c_t)^2)}
    } else if(OPTIMCONV == F){ # To ensure that OPTIMCONV exists - the user knows whether they want optimized convolutions or not.
      c_t = C_T(timeFrames[i+1],phiVec,thetaVec)
    }

    c_t = abs(c_t)
    if(DEBUG){cat("ct:",c_t,"at time",timeFrames[i+1],"\n")}
    #if(c_t<0){stop("DONESO")}
    logLHood = logLHood - (1/2)*exp(b)*((timeFrames[i+1]-timeFrames[i])/c_t)*((data[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    #logLHood = logLHood - (lambda/2)*((timeFrames[i+1]-timeFrames[i])/c_t)*((data[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    # if(DEBUG){cat("logLHood:",logLHood,"\n")}
    # if(DEBUG){cat('\n','Current LogLikelihood value ', i,'/', length(data) , ' : ',logLHood,'\n')}

  }

  logLHood = logLHood + length(data)*(b-log(2*pi))*(0.5)

  #debuging code
  # if(DEBUG){cat('\n','LogLikelihod calculated:',logLHood,'\n')}

  return(logLHood)
}


logLikelihoodLambda = function(data,timeFrames,phiVec,thetaVec,lambda){

  # We will need a summand as we are looking at the log transformed liklihood.
  logLHood = 0

  for(i in 1:length(data)){
    # We compute the tissue concentration(C_T - at the end point of each time frame)


    if(OPTIMCONV == T){
      c_t = intpolFitValue(i,phiVec,thetaVec)
      #if(c_t==0){return<--Inf}else{return<- return - (lambda/2)*((tVec[i+1]-tVec[i])/c_t)*((dataVec[i]-c_t)^2)}
    } else if(OPTIMCONV == F){ # To ensure that OPTIMCONV exists - the user knows whether they want optimized convolutions or not.
      c_t = C_T(timeFrames[i+1],phiVec,thetaVec)
    }

    c_t = abs(c_t)
    # if(DEBUG){cat("ct:",c_t,"at time",timeFrames[i+1],"\n")}
    #
    # resid = (data[i]-c_t)
    # if(DEBUG){cat("resid:",resid,"at time",timeFrames[i+1],"\n")}
    #
    #
    # adding = - (1/2)*lambda*((timeFrames[i+1]-timeFrames[i])/c_t)*((data[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    # if(DEBUG){cat("adding:",adding,"at time",timeFrames[i+1],"\n")}

    #if(c_t<0){stop("DONESO")}
    logLHood = logLHood - (1/2)*lambda*((timeFrames[i+1]-timeFrames[i])/c_t)*((data[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    #logLHood = logLHood - (lambda/2)*((timeFrames[i+1]-timeFrames[i])/c_t)*((data[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    # if(DEBUG){cat("logLHood:",logLHood,"\n")}
    # if(DEBUG){cat('\n','Current LogLikelihood value ', i,'/', length(data) , ' : ',logLHood,'\n')}

  }

  logLHood = logLHood + length(data)*(log(lambda/(2*pi)))*(0.5)
  # logLHood = logLHood + length(data)*(b-log(2*pi))*(0.5)

  #debuging code
  # if(DEBUG){cat('\n','LogLikelihod calculated:',logLHood,'\n')}

  return(logLHood)
}


petLogLikelihoodDensity = function(datum, parameters, dimParameters, timeFrames){
  modelOrder = (dimParameters-1)/2
    phiVec = parameters[1:modelOrder]
    thetaVec = parameters[(modelOrder+1):(2*modelOrder)]
    b = parameters[(2*modelOrder)+1]
  # We will need a summand as we are looking at the log transformed liklihood.
  logLHood = 0

  for(i in 1:length(datum)){
    # We compute the tissue concentration(C_T - at the end point of each time frame)


    c_t = intpolFitValue(i+1,phiVec,thetaVec)
    # if(OPTIMCONV == T){
    #   c_t = intpolFitValue(i+1,phiVec,thetaVec)
    #   #if(c_t==0){return<--Inf}else{return<- return - (lambda/2)*((tVec[i+1]-tVec[i])/c_t)*((dataVec[i]-c_t)^2)}
    # } else if(OPTIMCONV == F){ # To ensure that OPTIMCONV exists - the user knows whether they want optimized convolutions or not.
    #   c_t = C_T(timeFrames[i+1],phiVec,thetaVec)
    # }

    c_t = abs(c_t)
    # if(DEBUG){cat("ct:",c_t,"at time",timeFrames[i+1],"\n")}
    #if(c_t<0){stop("DONESO")}
    logLHood = logLHood - (1/2)*exp(b)*((timeFrames[i+1]-timeFrames[i])/c_t)*((datum[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    #logLHood = logLHood - (lambda/2)*((timeFrames[i+1]-timeFrames[i])/c_t)*((data[i]-c_t)^2) + (1/2)*log((timeFrames[i+1]-timeFrames[i])/c_t)
    # if(DEBUG){cat("logLHood:",logLHood,"\n")}
    # if(DEBUG){cat('\n','Current LogLikelihood value ', i,'/', length(datum) , ' : ',logLHood,'\n')}

  }

  logLHood = logLHood + length(datum)*(b-log(2*pi))*(0.5)

  #debuging code
  # if(DEBUG){cat('\n','LogLikelihod calculated:',logLHood,'\n')}

  return(logLHood)
}


petLogBioPriorDensity <- function(parameters, dimParameters,timeFrames){
  modelOrder = (dimParameters-1)/2
  phi = parameters[1:modelOrder]
  theta = parameters[(modelOrder+1):(2*modelOrder)]
  b = parameters[(2*modelOrder)+1]
  if(exp(b) == 0){
    return(log(0))
  }
    for(m in 1:modelOrder){
      if(phivPrior(phi[m]) == 0 | thetavPrior(theta[m]) == 0){
        # if(DEBUG){print('\n  NEGATIVE PROP VALUE \n')}
        # in this case the propValue is NOT in the prior so fails the propVPriorCheck
        return(log(0))
      }
    }
  # cat("log(dBioPrior(phi,theta)):",log(dBioPrior(phi,theta)))
  # cat("(alpha*b)-(beta*exp(b)):",(alpha*b)-(beta*exp(b)))
  (alpha*b)-(beta*exp(b)) + log(dBioPrior(phi,theta))
}

petLogVagueBioPriorDensity <- function(parameters, dimParameters,timeFrames){
  modelOrder = (dimParameters-1)/2
  phi = parameters[1:modelOrder]
  theta = parameters[(modelOrder+1):(2*modelOrder)]
  b = parameters[(2*modelOrder)+1]
  if(exp(b) == 0){
    return(log(0))
  }
  for(m in 1:modelOrder){
    if(phivPrior(phi[m]) == 0 | thetavPrior(theta[m]) == 0){
      # if(DEBUG){print('\n  NEGATIVE PROP VALUE \n')}
      # in this case the propValue is NOT in the prior so fails the propVPriorCheck
      return(log(0))
    }
  }
  # cat("log(dBioPrior(phi,theta)):",log(dBioPrior(phi,theta)))
  # cat("(alpha*b)-(beta*exp(b)):",(alpha*b)-(beta*exp(b)))
  (alpha*b)-(beta*exp(b))
}


logPosterior<-function(data,timeFrames,phi,theta,b){

  logPost = 0

  # if(DEBUG){cat('\n','Calculating Log-Likelihood....')}

  logLHood = logLikelihood(data,timeFrames,phi,theta,b) # the logLikeliHood function is show below

  # if(DEBUG){cat('\n','Multiplying Log-Likelihood by lambda prior...')}

  logPost = logLHood + (alpha*b)-(beta*exp(b))
  #logPost = logLHood + ((alpha - 1)*log(lambda)) - (beta*lambda)


  # if(DEBUG){cat('\n','Log-likelihod multiplied by lamda prior: ','\n',logPost)}


  #theta and phi priors - (we use vague priors here)
  # if(DEBUG){cat('\n ','Multiplying by theta/phi priors ...')}
  # for(i in 1:length(theta)){
  #   logPost = logPost + log(phivPrior(phi[i])*thetavPrior(theta[i]))
  #   if(DEBUG){cat('\n','Multiplying by', i , 'th theta/phi: ','\n',logPost)}
  # }

  logPost = logPost + log(dBioPrior(phi,theta))

  # if(DEBUG){cat('\n','Final LOG posterior output: ','\n',logPost,'\n')}


  return(logPost)

}


logPosteriorWAnn<-function(data,timeFrames,phi,theta,b,alphaT){

  logPost = 0

  # if(DEBUG){cat('\n','Calculating Log-Likelihood....')}

  logLHood = logLikelihood(data,timeFrames,phi,theta,b) # the logLikeliHood function is above

  # if(DEBUG){cat('\n','Multiplying Log-Likelihood by lambda prior...')}

  logPost = (logLHood*alphaT) + (alpha*b)-(beta*exp(b))
  #logPost = logLHood + ((alpha - 1)*log(lambda)) - (beta*lambda)


  # if(DEBUG){cat('\n','Log-likelihod multiplied by lamda prior: ','\n',logPost)}


  #theta and phi priors - (we use vague priors here)
  # if(DEBUG){cat('\n ','Multiplying by theta/phi priors ...')}
  # for(i in 1:length(theta)){
  #   logPost = logPost + log(phivPrior(phi[i])*thetavPrior(theta[i]))
  #   if(DEBUG){cat('\n','Multiplying by', i , 'th theta/phi: ','\n',logPost)}
  # }

  logPost = logPost + log(dBioPrior(phi,theta))

  # if(DEBUG){cat('\n','Final LOG posterior output: ','\n',logPost,'\n')}


  return(logPost)

}


petBioPriorSampler <- function(dimParameter, numSamples, timeFrames){
  particles = matrix(NA,numSamples,dimParameter)
  modelOrder = (dimParameter - 1)/2
  phiMean = 3e-3

  for(j in 1:modelOrder){
    particles[, j] = truncnorm::rtruncnorm(numSamples, a = 1e-5,b = 1e-2, mean = phiMean, sqrt(1e-3))

    if(j ==1){
      #samples[,1,j+modelId] = runif(numSamples,2e-4,1e-2)
      for(i in 1:numSamples){
        particles[i, (j+modelOrder)] = truncnorm::rtruncnorm(1, a = 2e-4, b = 1e-2,
                                                                  mean =particles[i,j]/15,
                                                                  sqrt(1e-2))
      }
      #samples[,1,(j+modelId)] = rtruncnorm(numSamples,a = 2e-4,b = 1e-2,mean =samples[,1,j]/15 ,1e-2)
      phiMean = 1e-3
    }else if(j == 2){
      #samples[,1,j+modelId] = runif(numSamples,samples[,1,j+modelId-1],6e-2)
      for(i in 1:numSamples){
        particles[i, (j+modelOrder)] = truncnorm::rtruncnorm(1, a = particles[i,(j+modelOrder-1)],
                                                                  b = 6e-2,
                                                                  mean =particles[i, j]/4,
                                                                  sqrt(1e-2))
      }
      # samples[,1,j+modelId] = rtruncnorm(numSamples,a = samples[,1,(j+modelId-1)],b = 6e-2,mean =samples[,1,j]/4 ,1e-2)
      cat("2 comp model")
    }else if( j == 3){
      #samples[,1,j+modelId] = runif(numSamples,samples[,1,j+modelId-1],6e-2)
      # stop("Fix the third intial value sampling!!")
      for(i in 1:numSamples){

        particles[i,(j+modelOrder)] = truncnorm::rtruncnorm(1,a = particles[i, (j+modelOrder-1)],b = 6e-2,mean =particles[i,j],sqrt(1e-2))

      }
      cat("3 comp model")
    }


  }

  particles[,2*modelOrder+1]= log(rgamma(numSamples,10^-3,rate = 10^-3)+1e-300)

  particles
}
