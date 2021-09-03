// [[Rcpp::depends(RcppArmadillo)]]
// #include <gsl/gsl_randist.h>
#include <iostream>
#include <string>
// #include <armadillo>
#include <cmath>
#include "truncated_normal.hpp"
#include "sourceVerOne.hpp"
#include <vector>
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace std::chrono;
using namespace Rcpp;

arma::mat rSamples;

// can make the paramters this into pointer.
arma::vec mcmcMove(arma::vec initialValue, int steps, int modelOrder, double alphaT, arma::vec data,arma::mat cholDecomp, double &acceptanceRate, arma::vec timeFrames){
  int dim  = initialValue.n_elem;

  arma::mat markovChain(dim,steps+1);
  markovChain.col(0) = initialValue;
  arma::vec proposalState(dim);

  // double acceptanceRate = 0;

  for(int i  = 1; i <steps+1;i++){
    // std::cout<< "\n Currently at step:"<< i <<" out of "<< steps<< " steps \n";
    // cout<< "\n markovChain.chol(i-1): \n"<<markovChain.col(i-1);
    // cout<< "\n chol times randn \n "<< arma::trans(cholDecomp)*arma::randn(dim);
    NumericVector innovations = rnorm(dim);

    std::vector<double> tempA = as<std::vector<double> >(innovations);
    arma::vec tempB = arma::conv_to<arma::vec>::from(tempA);

    // proposalState = markovChain.col(i-1) + arma::trans(cholDecomp)*arma::randn(dim);
    proposalState = markovChain.col(i-1) + arma::trans(cholDecomp)*tempB;
    // cout<< "\n proposalState \n"<< proposalState;
  // for efficiency we check if the proposed value is in the vague prior (i.e. within correct limit)
  double propLogPriorDensity = petLogBioPriorDensity(proposalState,dim);
   bool proposalVaguePriorCheck  = true;
   if(exp(propLogPriorDensity) ==0){
     // std::cout<< "\n Rejected:"<<proposalState;
     proposalVaguePriorCheck  = false;
   }

   // if(exp(proposalState(dim-1)) == 0 ){
   //   proposalVaguePriorCheck = false;}
   //
   // if(proposalVaguePriorCheck){
   //    for(int m = 0; m<modelOrder;m++){
   //      if(vaguePriorPhi(proposalState(m))==0 ||  vaguePriorTheta(proposalState(m+modelOrder)) ==0){
   //        proposalVaguePriorCheck  = false;
   //      }
   //    }
   //  }

   // cout<< "prior Check Complete, priorCheck is:"<<proposalVaguePriorCheck << "\n";
   double acceptanceRatio;

   if( proposalVaguePriorCheck == false){ acceptanceRatio = 0;}else{
     double proposalStateLogProbabiliy = logPosteriorWithAnnealing(data,timeFrames,
                                                                   proposalState.subvec(0,modelOrder-1),proposalState.subvec(modelOrder,(2*modelOrder)-1),proposalState(dim-1),alphaT);

      // cout<< "proposal state Log Probability calculated";
     double currentStateLogProbabiliy = logPosteriorWithAnnealing(data,timeFrames,
                                                                  markovChain.col(i-1).subvec(0,modelOrder-1),markovChain.col(i-1).subvec(modelOrder,(2*modelOrder)-1),markovChain.col(i-1)(dim-1),alphaT);

      // cout<< "current state log probability calculated";
   // Some state checks missing i.e. if YLogProb == -Inf etc etc...
     acceptanceRatio = exp(proposalStateLogProbabiliy-currentStateLogProbabiliy);
     // cout<< "\n Acceptance Ratio:"<< acceptanceRatio;
   }
   // cout<< "\n acceptance Ratio Calculated";

  // arma::vec u = arma::randu(1);
  NumericVector u = runif(1);
  if(u(0)<acceptanceRatio){
    // cout<< "\n Accepted";
    markovChain.col(i) = proposalState;
    acceptanceRate++;
  }else{
   // cout<< "\n Rejected";
    markovChain.col(i) = markovChain.col(i-1);
  }

  }

  // std::cout<<"\n MCMC Kerenl compeleted";

  acceptanceRate = acceptanceRate/steps;
  // stop("STOP");
  return(markovChain.col(steps));
}


arma::mat adaptiveStepSize(arma::mat samples){
  int dim = samples.n_cols;
  // cout<< "\n DIMMMMMMMMMMMMMM:"<<dim;
  arma::mat sigma = arma::cov(samples);
  // cout<< "\n Sigma: \n"<< sigma;
  arma::mat sigmaScaled = ((pow(2.3,2))/dim)*sigma;
  // arma::mat sigmaScaled(dim,dim);
  // sigmaScaled.fill(0);
  // sigmaScaled.diag().fill(1e-10);
  // if(arma::det(sigmaScaled)<10e-60){
  //   arma::mat add(dim,dim);
  //   add.fill(10e-40);
  //   sigmaScaled = sigmaScaled + add;
  // }
  return sigmaScaled;
}

arma::vec meanWeight(arma::mat samples, arma::vec weights){
  // cout<< "\n weights:"<< weights;
  // cout<< "\n SUM OF WEIGHTS : "<< arma::sum(weights);
  arma::vec ret(samples.n_cols);
  for(int i = 0; i<samples.n_cols;i++){
    double meanAti_j =0;
    for(int j = 0; j<samples.n_rows;j++){

      meanAti_j  = meanAti_j + (weights(j)*samples(j,i));
    }
    // cout<< "\n meanAti_j ;\n"<< meanAti_j;
    ret(i) = meanAti_j;
  }
  return(ret);
}

arma::mat covWeight(arma::mat samples, arma::vec weights){
  // exp(weights);
  arma::mat ret(samples.n_cols,samples.n_cols);
  ret.fill(0);
  arma::vec weightedMeans = meanWeight(samples,weights);

  for(int i =0; i<samples.n_cols;i++){
    for(int j = i;j<samples.n_cols;j++){
      double covarainceAti_j = 0;
      for(int k = 0;k<samples.n_rows;k++){
        covarainceAti_j = covarainceAti_j + weights(k)*((samples(k,i) - weightedMeans(i))*(samples(k,j) - weightedMeans(j)));
      }
      arma::vec weightsSquared = arma::pow(weights,2);
      ret(i,j) = covarainceAti_j/(1-arma::sum(weightsSquared));
    }
  }
  // cout<< "\n COmputed size:"<< ret.size();

  // cout<< "\n Ret before trans \n "<< ret;
  arma::mat transRet = arma::trans(ret);

  // cout<< "\n Ret after trans \n "<< ret;
  // cout<< "\n transRet after trans \n "<< transRet;
  // cout<< "\n diagmat(ret) \n"<<arma::diagmat(ret);
  ret = ret+ transRet - arma::diagmat(ret);

  // cout<< "\n Ret after mirror \n "<< ret;
  return(ret);

}

// not done correctly
arma::mat weightedAdaptiveStepSize(arma::mat samples, arma::vec logWeights){
  int dim = samples.n_cols;
  // arma::mat sigma = arma::var(samples,0,1);
  arma::mat sigma = covWeight(samples,exp(logWeights));
  arma::mat sigmaScaled = ((pow(2.3,2))/dim)*sigma;
  // std::cout<<"\n Sigma: \n"<< sigma;
  // std::cout<<"\n Sigma size: \n"<< sigma.size();
  // std::cout<< "\n Sigma scaled: \n"<<sigmaScaled;
  // std::cout<< "\n Sigma sclaed size: \n"<<sigmaScaled.size();
  return sigmaScaled;
}

// Container for full SMC particle+weights system
class SMCParticlesSystem{
  public:
    int dim,numOfSamples,numOfDistributions;
    // container for samples
    arma::cube samples;


    // all types of weights
    arma::mat logWeights;
    arma::mat normLogWeights;
    arma::mat unNormLogWeights;
    arma::mat incremWeights;
    arma::mat barLogWeights;
    arma::umat reSamplingIndices;

    //normalising constant
    double normalisingConstant;

    //other metadata/quantities
    arma::vec ESS;
    vector<arma::mat> optimialStepSizes;
    arma::vec alphaVec;
    arma::vec accRate;


    SMCParticlesSystem(int newNumOfSamples,int newNumOfDistributions,int newDim){
      // cout<< "\n creating particles";
      dim = newDim;
      numOfSamples = newNumOfSamples;
      numOfDistributions = newNumOfDistributions;
      samples.resize(numOfSamples,numOfDistributions,dim);

      logWeights.resize(numOfSamples,numOfDistributions);
      normLogWeights.resize(numOfSamples,numOfDistributions);
      unNormLogWeights.resize(numOfSamples,numOfDistributions);
      incremWeights.resize(numOfSamples,numOfDistributions);
      barLogWeights.resize(numOfSamples,numOfDistributions);
      reSamplingIndices.resize(numOfSamples,1);

      alphaVec.resize(numOfDistributions);
      ESS.resize(numOfDistributions);
      accRate.resize(numOfDistributions);
    }
    void initialiseParticlesWithPrior();
    void initialiseWeights();
    void reSample(int currentDistribution);
    void initialiseParticlesWithR();
    void addResamples(arma::uvec reSampledIndices);
};

void SMCParticlesSystem::addResamples(arma::uvec reSampledIndices){
  // cout<<"\m reSamplingIndices.n_cols" ;
  reSamplingIndices.insert_cols((reSamplingIndices.n_cols),1);
  reSamplingIndices.col(reSamplingIndices.n_cols-1) = reSampledIndices;

}
void SMCParticlesSystem::initialiseParticlesWithPrior(){
  // cout<< "\n Initialising particles \n ";
  // cout<< "\n dim: "<< dim;
  int modelOrder = (dim-1)/2;
  // cout<< "\n Model Order: "<< modelOrder;

  for(int j = 1; j<=modelOrder;j++){
    // cout<< "\n j: "<<j;
    for(int i = 0;i<samples.n_rows;i++){

      arma::vec u = arma::randu(1);
      int seed = abs(u(0))*100 + 1;

      // double truncated_normal_ab_sample ( double mu, double sigma, double a, double b, int &seed );
      // samples(i,j) = truncated_normal_ab_sample(3e-3,sqrt(1e-3),1e-5,1e-2,seed);

      if(j==1){
        samples(i,0,j-1) = truncated_normal_ab_sample(3e-3,sqrt(1e-3),1e-5,1e-2,seed);
        samples(i,0,j+modelOrder-1) = truncated_normal_ab_sample(samples(i,0,j)/15,sqrt(1e-2),2e-4,1e-2,seed);
        // cout<< "samples(i,0,j)"<<samples(i,0,j);
      }else if(j==2)
      {
        samples(i,0,j-1) = truncated_normal_ab_sample(1e-3,sqrt(1e-3),1e-5,1e-2,seed);
        samples(i,0,j+modelOrder-1) = truncated_normal_ab_sample(samples(i,0,j)/4,sqrt(1e-2),samples(i,0,j+modelOrder-1),6e-2,seed);

      } else if(j==3)
      {
        samples(i,0,j-1) = truncated_normal_ab_sample(1e-3,sqrt(1e-3),1e-5,1e-2,seed);
        samples(i,0,j+modelOrder-1) = truncated_normal_ab_sample(samples(i,0,j),sqrt(1e-2),samples(i,0,j+modelOrder-2),6e-2,seed);
      }
    }
    samples(arma::span::all,arma::span(0),arma::span(dim-1)) = log(arma::randg<arma::vec>(samples.n_rows,arma::distr_param(10e-3,1/10e-3)) +1e-300) ;
    // samples(arma::span::all,arma::span(0),arma::span(dim-1)) = log(rgamma(numOfSamples,10^-3,10^-3)+1e-300) ;
    // samples(arma::span::all,arma::span(0),arma::span(dim-1)).fill(-250);
      }
}

void SMCParticlesSystem::initialiseWeights(){
  incremWeights(arma::span::all,arma::span(0,0)).fill(log(1));
  logWeights(arma::span::all,arma::span(0,0)).fill(-log(numOfSamples));
}


// The SMC sampler
SMCParticlesSystem smc(int numOfSamples, int numOfDistrubtions, arma::vec data, arma::vec timeFrames, int modelOrder){

  SMCParticlesSystem particleCollection(numOfSamples,numOfDistrubtions,(2*modelOrder)+1);

  // std::cout<< "\n Particles, weights initialising";
  // particleCollection.initialiseParticlesWithPrior();
  particleCollection.initialiseParticlesWithR();

  particleCollection.initialiseWeights();
  // cout<< "\n -log(num)"<< -log(numOfSamples);
  // cout<< "\n sum of logweights:"<<arma::sum(exp(particleCollection.logWeights.col(0)));

  // std::cout<< "\n Particles, weights initialised";

  particleCollection.alphaVec(0) = 0;

  particleCollection.accRate(0) =0;

  arma::uvec numSamplesIndices = arma::linspace<arma::uvec>(0,numOfSamples-1,numOfSamples);

  arma::vec normWeightsAtZero = exp(particleCollection.logWeights.col(0));
  // cout<< "\n ESS WEGHT CALUC SUM:"<< arma::sum(normWeights);
  arma::vec normWeightsSquaredAtZero  = arma::pow(normWeightsAtZero,2);
  particleCollection.ESS(0) = 1.0/arma::sum(normWeightsSquaredAtZero);

  for(int interimDistIndex = 1; interimDistIndex<numOfDistrubtions;interimDistIndex++){
    // cout<<"\n Currently at step: "<< interimDistIndex+1;
    double binSearchOutput = pow((double(interimDistIndex+1)/numOfDistrubtions),5);
    // std::cout<< "\n bin search output:"<<binSearchOutput;
    // std::cout<< "\n pow(1/500,5)"<<pow(1.0/500,5);
    // std::cout<< "\n 1/500"<<1.0/50;
    // if(interimDistIndex == 100){break;}
    // break;
    particleCollection.alphaVec(interimDistIndex) = binSearchOutput;
    bool resamplingCheck  = false;
    arma::mat stepSize;

    // double ESS = 1/arma::sum(pow(exp(particleCollection.logWeights.col(interimDistIndex-1)),2));
    // double ESS = 10;
    // std::cout<< "\n ESS:"<< ESS;
    // cout<< "\n Log weights:\n "<<particleCollection.logWeights.col(interimDistIndex-1);
    arma::vec normWeights = exp(particleCollection.logWeights.col(interimDistIndex-1));
    // cout<< "\n ESS WEGHT CALUC SUM:"<< arma::sum(normWeights);
    arma::vec normWeightsSquared  = arma::pow(normWeights,2);
    particleCollection.ESS(interimDistIndex) = 1.0/arma::sum(normWeightsSquared);
    // if(interimDistIndex == 3){break;}
    // particleCollection.ESS(interimDistIndex-1) = 1.0/arma::sum(pow(exp(particleCollection.logWeights.col(interimDistIndex-1)),2));
    // particleCollection.ESS(interimDistIndex-1) = ESS;

    // std::cout<< "\n ESS Calculated: "<< particleCollection.ESS(interimDistIndex-1);
    // bool check = particleCollection.ESS(interimDistIndex-1)<0.5*numOfSamples;
    // // bool check = true;
    // std::cout<< "\n Check: "<< check;

    if(particleCollection.ESS(interimDistIndex)<0.5*numOfSamples){
      //Resample if needed
      // std::cout<< "\n Resampling";
      // std::cout<< "\n numSamplesIndices size"<< numSamplesIndices.size();
      // std::cout<< "\n particleCollection.logWeights.col(interimDistIndex-1) size:"<< particleCollection.logWeights.n_rows;
      // arma::uvec reSampledIndices = sample_main(numSamplesIndices,numOfSamples,true,exp(particleCollection.logWeights.col(interimDistIndex-1)));
      arma::uvec reSampledIndices =  csample(numOfSamples, exp(particleCollection.logWeights.col(interimDistIndex-1)));
      particleCollection.addResamples(reSampledIndices);
      // arma::uvec reSampledIndices ;
      // cout<< "resampled nices : \n"<<reSampledIndices;
      arma::mat previousParticles = particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all);
      // cout<< "PREV PARTICLES:"<<previousParticles(arma::span(1),arma::span::all);
      // cout<< "PREV PARTICLES ROWS : \n"<<previousParticles.rows(reSampledIndices);
      particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all) = previousParticles.rows(reSampledIndices);
      // cout<< "particle collection out \n: " << particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all);
      particleCollection.barLogWeights.col(interimDistIndex-1).fill(-log(numOfSamples));
      stepSize = adaptiveStepSize(particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all));
      // previousParticles = particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all);
      // stepSize = weightedAdaptiveStepSize(previousParticles,exp(particleCollection.logWeights.col(interimDistIndex-1)));
      // break;
    }else{
      // std::cout<< "\n No resampling";
      particleCollection.barLogWeights.col(interimDistIndex-1) = particleCollection.logWeights.col(interimDistIndex-1);
      arma::mat previousParticles = particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all);
      // cout<< "\n previousParticles.n_cols: "<<previousParticles.n_cols;
      // cout<< "\n previousParticles.n_rows: "<<previousParticles.n_rows;
      // break;
      // stepSize = adaptiveStepSize(particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all));

      // cout<<"\n logWeights:"<< particleCollection.logWeights.col(interimDistIndex-1);
      // cout<<"\n exp(logWeights):"<< exp(particleCollection.logWeights.col(interimDistIndex-1));
      // cout<<"\n logWeights:"<< particleCollection.logWeights.col(interimDistIndex-1);
      stepSize = weightedAdaptiveStepSize(previousParticles,particleCollection.logWeights.col(interimDistIndex-1));
      // break;
      // stepSize = adaptiveStepSize(particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all));
    }
    // particleCollection.optimialStepSizes[interimDistIndex] = stepSize;
    // std::cout<< " \n New step size computed: \n"<< stepSize;
    // break;
    double det = arma::det(stepSize);
    // std::cout<< " \n Determinant of stepsize:"<< det;
    arma::mat stepSizeCholDecomp = arma::chol(stepSize);
    // std::cout<< "\n Chol Decomp StepSize computed"<<stepSizeCholDecomp;

    for(int particleIndex =0; particleIndex< numOfSamples;particleIndex++){

      arma::vec previousPhi = particleCollection.samples(arma::span(particleIndex,particleIndex),arma::span(interimDistIndex-1,interimDistIndex-1),arma::span(0,modelOrder-1));
      arma::vec previousTheta = particleCollection.samples(arma::span(particleIndex,particleIndex),arma::span(interimDistIndex-1,interimDistIndex-1),arma::span(modelOrder,2*modelOrder-1));
      double previousB = particleCollection.samples(particleIndex,interimDistIndex-1,2*modelOrder);

      double previousSampleLogLikelihood = logLikelihood(data,timeFrames,previousPhi,previousTheta,previousB);

      // std::cout<< "\n Previous Sample Log likelihood computed:"<<previousSampleLogLikelihood;
      // std::cout<< "\n particleCollection.alphaVec(interimDistIndex)"<<particleCollection.alphaVec(interimDistIndex);
      // std::cout<< "\n particleCollection.alphaVec(interimDistIndex-1)"<<particleCollection.alphaVec(interimDistIndex-1);

      particleCollection.incremWeights(particleIndex,interimDistIndex) = (particleCollection.alphaVec(interimDistIndex)*previousSampleLogLikelihood) - (particleCollection.alphaVec(interimDistIndex-1)*previousSampleLogLikelihood);
      // std::cout<< "\n Incremental weights:"<<particleCollection.incremWeights(particleIndex,interimDistIndex);
      particleCollection.logWeights(particleIndex,interimDistIndex) = particleCollection.barLogWeights(particleIndex,interimDistIndex-1)+particleCollection.incremWeights(particleIndex,interimDistIndex);

      double particleAcceptanceRate =0;
      arma::vec propogatedParticles = mcmcMove(particleCollection.samples(arma::span(particleIndex,particleIndex),arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all)
                                                 ,3,modelOrder,particleCollection.alphaVec(interimDistIndex),data,stepSizeCholDecomp,particleAcceptanceRate, timeFrames);
      particleCollection.accRate(interimDistIndex) = particleCollection.accRate(interimDistIndex)+particleAcceptanceRate;
      // std::cout<< "\n MCMC Kernel Completed";
      particleCollection.samples(arma::span(particleIndex,particleIndex),arma::span(interimDistIndex,interimDistIndex),arma::span::all) = propogatedParticles;
    }
    // std::cout<< "\n All Particles propogated";
    particleCollection.accRate(interimDistIndex) =  particleCollection.accRate(interimDistIndex)/numOfSamples;
    // std::cout<< "\n Acceptance rate:"<<particleCollection.accRate(interimDistIndex) ;

    particleCollection.unNormLogWeights.col(interimDistIndex) = particleCollection.logWeights.col(interimDistIndex);
    double maxLogWeight = arma::max(particleCollection.logWeights.col(interimDistIndex));
    particleCollection.logWeights.col(interimDistIndex) = particleCollection.logWeights.col(interimDistIndex) - maxLogWeight - log(arma::sum(exp(particleCollection.logWeights.col(interimDistIndex)-maxLogWeight)));
  }

  // cout<< "\n SMC Completed";
  return particleCollection;
}

//'Cpp SMC Sampler Implementation (PET ONLY)
//'
//'
//'IMPORTANT: Recommended that you use \code{smcSamplerRcpp} instead !!
//'IMPORTANT: This function is used for development and unit testing only!!
//'\code{smcSamplerCpp} Returns a List containing various values of the SMC
//'sampler e.g.  particle and (log) weight collection, normalising constant of
//'posterior denisty and other adaptive techniques metadata.
//'
//'This is NOT a general function  like the other samplers, instead this works only
//'for the PET model. However, like \code{smcSampler}, it will provide weighted
//'samples that can be used to approximate different properties of the Bayesian
//'model of interest
//'
//'@param numParticles A number specifiying the number of particles to be used in
//'  the SMC sampler.
//'@param numOfDistributions A number specifiying the number of intermediate
//'  distributions to be used in the simulated annealing (when adpatice CESS is
//'  not used).
//'@param dataset A vector containing the data point.
//'@param timeFrames A vector containing the data point.
//'@param dimParameter The dimension of the paramaters space.
//'
//'@export
// [[Rcpp::export]]
List smcSamplerCpp(int numParticles,int numOfDistrubtions, arma::vec dataSet, arma::vec timeFrames, int dimParameters){
  int modelOrder = (dimParameters-1)/2;
  auto start = high_resolution_clock::now();
  SMCParticlesSystem  testingSMC = smc(numParticles,numOfDistrubtions,dataSet,timeFrames,modelOrder);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout << duration.count() << endl;
  double normalizingConstant = arma::mean(exp(testingSMC.incremWeights.col(0)));
    for(int i = 1; i<numOfDistrubtions;i++){
      normalizingConstant = normalizingConstant * arma::sum(exp(testingSMC.incremWeights.col(i)+testingSMC.barLogWeights.col(i-1)));
    }
    // return(normalizingConstant);
  return(List::create(Named("particles") = testingSMC.samples, Named("normLogWeights") = testingSMC.logWeights, Named("UnNormlogWeights") = testingSMC.unNormLogWeights,
  Named("alphaVec") = testingSMC.alphaVec, Named("acceptanceRates") = testingSMC.accRate, Named("incremWeights") = testingSMC.incremWeights,
  Named("barLogWeights") = testingSMC.barLogWeights, Named("ESS") = testingSMC.ESS, Named("normalisingConstant") = normalizingConstant,
  Named("reSamplingIndex") = testingSMC.reSamplingIndices));

           // // container for samples
           // arma::cube samples;
           //
           //
           // // all types of weights
           // arma::mat logWeights;
           // arma::mat normLogWeights;
           // arma::mat unNormLogWeights;
           // arma::mat incremWeights;
           // arma::mat barLogWeights;
           //
           // //normalising constant
           // double normalisingConstant;
           //
           // //other metadata/quantities
           // arma::vec ESS;
           // vector<arma::mat> optimialStepSizes;
           // arma::vec alphaVec;
           // arma::vec accRate;
}

// int main(){
  // int arg =0 ;
  // double test = smcSamplerRcpp(200,500, data,timeFrames, 5);
  // cout<<test;


  // return(arg);
// }

// [[Rcpp::export]]
void rTakeit(arma::mat test){
  rSamples  = test;
  // std::cout<< "\n rSampels:"<<rSamples;
}

void SMCParticlesSystem::initialiseParticlesWithR(){
  samples(arma::span::all,arma::span(0),arma::span::all) = rSamples;
}

// [[RRcpp::export]]
// arma::vec sampleRnorm(){
  // NumericVector ret = rnorm(20);
  // std::vector<double> b = as<std::vector<double> >(ret);
  // arma::vec A = arma::conv_to<arma::vec>::from(b);
  // std::cout<< "\n rSampels:"<<ret;
  // std::cout<< "\n rSampels:"<<b[11];
  // std::cout<< "\n rSampels:"<<A;
  // return(A);
// }
