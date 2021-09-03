// [[Rcpp::depends(RcppArmadillo)]]
// #include <gsl/gsl_randist.h>
#include <iostream>
#include <string>
// #include <armadillo>
#include <cmath>
// #include "truncated_normal.hpp"
#include "sourceVerOne.hpp"
#include <vector>
#include <RcppArmadillo.h>
#include <Rmath.h>
// using namespace std::chrono;
using namespace Rcpp;
using namespace std;

arma::vec rToySamples;

double toyLogLikelihood(double data, double parameters,double likeSigma){
  return R::dnorm(data, parameters, likeSigma, true);
}
// this doesntWork for some reason???
double toyLogBioPriorDensity(double parameters, double mu, double priorSigma){
  return R::dnorm(parameters, mu, priorSigma, true);
}

// [[Rcpp::export]]
double mcmcMoveToy(double initialValue,
                   int steps,
                   double alphaT,
                   double data,
                   double stepSize,
                   double priorSigma,
                   double mu,
                   double likeSigma){

                   // double &acceptanceRate,
  arma::vec markovChain(steps+1);
  markovChain(0) = initialValue;
  double proposalState;

  // double acceptanceRate = 0;

  for(int i  = 1; i <steps+1;i++){
    // std::cout<< "\n Currently at step:"<< i <<" out of "<< steps<< " steps \n";
    // cout<< "\n markovChain.chol(i-1): \n"<<markovChain.col(i-1);
    // cout<< "\n chol times randn \n "<< arma::trans(cholDecomp)*arma::randn(dim);
    double innovations = rnorm(1,0,1)[0]*stepSize;

    // proposalState = markovChain.col(i-1) + arma::trans(cholDecomp)*arma::randn(dim);
    proposalState = markovChain(i-1) + innovations;
    // cout<< "\n proposalState \n"<< proposalState;
  // for efficiency we check if the proposed value is in the vague prior (i.e. within correct limit)
  double propLogPriorDensity = R::dnorm(proposalState, mu, priorSigma, true);
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

     double currLogLikelihood = toyLogLikelihood(data, markovChain(i-1), likeSigma);
     double propLogLikelihood = toyLogLikelihood(data, proposalState, likeSigma);
     // double currLogPriorDensity = toyLogBioPriorDensity(markovChain(i-1), mu, priorSigma);
     double currLogPriorDensity = R::dnorm(markovChain(i-1), mu, priorSigma, true);

     // cout<< "currLogPriorDensity:"<<currLogPriorDensity;

     double proposalStateLogProbabiliy =  (propLogLikelihood * alphaT) + propLogPriorDensity;

      // cout<< "proposal state Log Probability calculated";
     double currentStateLogProbabiliy = (currLogLikelihood * alphaT) + currLogPriorDensity;

     // markovChain(i-1)double currentStateLogProbabiliy = `logPosteriorWithAnnealing(data,timeFrames,
                                                                  // markovChain.col(i-1).subvec(0,modelOrder-1),markovChain.col(i-1).subvec(modelOrder,(2*modelOrder)-1),markovChain.col(i-1)(dim-1),alphaT);`
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
    markovChain(i) = proposalState;
    // acceptanceRate++;
  }else{
   // cout<< "\n Rejected";
    markovChain(i) = markovChain(i-1);
  }

  }

  // std::cout<<"\n MCMC Kerenl compeleted";

  // acceptanceRate = acceptanceRate/steps;
  // cout<< " \n markovChain(i):"<<markovChain(steps);
  // stop("Stopped");// stop("STOP");
  return(markovChain(steps));
}


double adaptiveStepSizeToy(arma::vec samples){
  //int dim = samples.n_cols;
  // cout<< "\n DIMMMMMMMMMMMMMM:"<<dim;
  arma::mat sigma = arma::cov(samples);
  // cout<< "\n Sigma: \n"<< sigma;
  double sigmaScaled = ((pow(2.3,2))/1)*sigma(0,0);
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

// double meanWeightToy(arma::vec samples, arma::vec weights){
//   // cout<< "\n weights:"<< weights;
//   // cout<< "\n SUM OF WEIGHTS : "<< arma::sum(weights);
//   double ret;
//   for(int i = 0; i<samples.;i++){
//     double meanAti_j =0;
//     for(int j = 0; j<samples.n_rows;j++){
//
//       meanAti_j  = meanAti_j + (weights(j)*samples(j,i));
//     }
//     // cout<< "\n meanAti_j ;\n"<< meanAti_j;
//     ret(i) = meanAti_j;
//   }
//   return(ret);
// }

double covWeightToy(arma::vec samples, arma::vec weights){
  // exp(weights);
  double ret =0 ;
  double weightedMeans = dot(samples,weights)/arma::sum(weights);

  for(int i =0; i<samples.n_elem;i++){
    ret += weights(i)*(pow(samples(i)-weightedMeans,2));
  }

  return(ret/arma::sum(weights));

}

// not done correctly
double weightedAdaptiveStepSizeToy(arma::vec samples, arma::vec logWeights){
  //int dim = samples.n_cols;
  // arma::mat sigma = arma::var(samples,0,1);
  double sigma = covWeightToy(samples,exp(logWeights));
  double sigmaScaled = ((pow(2.3,2))/1)*sigma;
  // std::cout<<"\n Sigma: \n"<< sigma;
  // std::cout<<"\n Sigma size: \n"<< sigma.size();
  // std::cout<< "\n Sigma scaled: \n"<<sigmaScaled;
  // std::cout<< "\n Sigma sclaed size: \n"<<sigmaScaled.size();
  return sigmaScaled;
}

// Container for full SMC particle+weights system
class SMCToyParticlesSystem{
  public:
    int dim,numOfSamples,numOfDistributions;
    // container for samples
    arma::mat samples;


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


    SMCToyParticlesSystem(int newNumOfSamples,int newNumOfDistributions,int newDim){
      // cout<< "\n creating particles";
      dim = newDim;
      numOfSamples = newNumOfSamples;
      numOfDistributions = newNumOfDistributions;
      samples.resize(numOfSamples,numOfDistributions);

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
    // void initialiseParticlesWithR();
    void initialiseParticlesWithRToy();
    void addResamples(arma::uvec reSampledIndices);
};

void SMCToyParticlesSystem::addResamples(arma::uvec reSampledIndices){
  // cout<<"\m reSamplingIndices.n_cols" ;
  reSamplingIndices.insert_cols((reSamplingIndices.n_cols),1);
  reSamplingIndices.col(reSamplingIndices.n_cols-1) = reSampledIndices;

}

void SMCToyParticlesSystem::initialiseWeights(){
  incremWeights(arma::span::all,arma::span(0,0)).fill(log(1));
  logWeights(arma::span::all,arma::span(0,0)).fill(-log(numOfSamples));
}



void SMCToyParticlesSystem::initialiseParticlesWithRToy(){
  // cout<<"\n rToySamples:"<<rToySamples;
  samples(arma::span::all,0) = rToySamples;
}



// The SMC sampler

SMCToyParticlesSystem smcToy(int numOfSamples, int numOfDistrubtions, double data,
                          int dimParameters, double priorSigma, double likeSigma,
                          double mu){

  SMCToyParticlesSystem particleCollection(numOfSamples,numOfDistrubtions,dimParameters);

  // std::cout<< "\n Particles, weights initialising";
  // particleCollection.initialiseParticlesWithPrior();
  particleCollection.initialiseParticlesWithRToy();
  // std::cout<< "\n Prior Sample Particles:"<< particleCollection.samples(arma::span::all,0) ;
  particleCollection.initialiseWeights();
  // cout<< "\n -log(num)"<< -log(numOfSamples);
  // cout<< "\n sum of logweights:"<<arma::sum(exp(particleCollection.logWeights.col(0)));

  // std::cout<< "\n Particles, weights initialised";

  particleCollection.alphaVec(0) = 0;

  particleCollection.accRate(0) =0;

  arma::uvec numSamplesIndices = arma::linspace<arma::uvec>(0,numOfSamples-1,numOfSamples);

  arma::vec normWeightsAtZero = exp(particleCollection.logWeights.col(0));
  arma::vec normWeightsSquaredAtZero  = arma::pow(normWeightsAtZero,2);
  particleCollection.ESS(0) = 1.0/arma::sum(normWeightsSquaredAtZero);
  // stop("ISTOPED IT NOW");
  // abort();
  // cout<< "\n ESS WEGHT CALUC SUM:"<<particleCollection.ESS(0);
  for(int interimDistIndex = 1; interimDistIndex<numOfDistrubtions;interimDistIndex++){
    // cout<<"\n Currently at step: "<< interimDistIndex+1;
    double binSearchOutput = pow((double(interimDistIndex+1)/numOfDistrubtions),5);
    // std::cout<< "\n bin search output:"<<binSearchOutput;
    // std::cout<< "\n pow(1/500,5)"<<pow(1.0/500,5);
    // std::cout<< "\n 1/500"<<1.0/50;
    // if(interimDistIndex == 100){break;}
    // break;
    particleCollection.alphaVec(interimDistIndex) = binSearchOutput;
    // bool resamplingCheck  = false;
    double stepSize;

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

    // Adaptive Resampling
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
      arma::mat previousParticles = particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1));
      // cout<< "PREV PARTICLES:"<<previousParticles(arma::span(1),arma::span::all);
      // cout<< "PREV PARTICLES ROWS : \n"<<previousParticles.rows(reSampledIndices);
      particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1)) = previousParticles.rows(reSampledIndices);
      // cout<< "particle collection out \n: " << particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all);
      particleCollection.barLogWeights.col(interimDistIndex-1).fill(-log(numOfSamples));
      stepSize = adaptiveStepSizeToy(particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1)));
      // previousParticles = particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all);
      // stepSize = weightedAdaptiveStepSize(previousParticles,exp(particleCollection.logWeights.col(interimDistIndex-1)));
      // break;
    }else{
      // std::cout<< "\n No resampling";
      particleCollection.barLogWeights.col(interimDistIndex-1) = particleCollection.logWeights.col(interimDistIndex-1);
      arma::mat previousParticles = particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1));
      // cout<< "\n previousParticles.n_cols: "<<previousParticles.n_cols;
      // cout<< "\n previousParticles.n_rows: "<<previousParticles.n_rows;
      // break;
      // stepSize = adaptiveStepSize(particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all));

      // cout<<"\n logWeights:"<< particleCollection.logWeights.col(interimDistIndex-1);
      // cout<<"\n exp(logWeights):"<< exp(particleCollection.logWeights.col(interimDistIndex-1));
      // cout<<"\n logWeights:"<< particleCollection.logWeights.col(interimDistIndex-1);
      stepSize = weightedAdaptiveStepSizeToy(previousParticles,particleCollection.logWeights.col(interimDistIndex-1));
      // break;
      // stepSize = adaptiveStepSize(particleCollection.samples(arma::span::all,arma::span(interimDistIndex-1,interimDistIndex-1),arma::span::all));
    }

    // cout<<"stepSize"<<stepSize;
    // stop("stoped");

    // particleCollection.optimialStepSizes[interimDistIndex] = stepSize;
    // std::cout<< " \n New step size computed: \n"<< stepSize;
    // break;
    // double det = arma::det(stepSize);
    // std::cout<< " \n Determinant of stepsize:"<< det;
    // arma::mat stepSizeCholDecomp = arma::chol(stepSize);
    // std::cout<< "\n Chol Decomp StepSize computed"<<stepSizeCholDecomp;

    for(int particleIndex =0; particleIndex< numOfSamples;particleIndex++){

      double previousSampleLogLikelihood = toyLogLikelihood(data,
                                                            particleCollection.samples(particleIndex,interimDistIndex-1),
                                                            likeSigma);

      // std::cout<< "\n Previous Sample Log likelihood computed:"<<previousSampleLogLikelihood;
      // std::cout<< "\n particleCollection.alphaVec(interimDistIndex)"<<particleCollection.alphaVec(interimDistIndex);
      // std::cout<< "\n particleCollection.alphaVec(interimDistIndex-1)"<<particleCollection.alphaVec(interimDistIndex-1);

      particleCollection.incremWeights(particleIndex,interimDistIndex) = (particleCollection.alphaVec(interimDistIndex)*previousSampleLogLikelihood) - (particleCollection.alphaVec(interimDistIndex-1)*previousSampleLogLikelihood);
      // std::cout<< "\n Incremental weights:"<<particleCollection.incremWeights(particleIndex,interimDistIndex);
      particleCollection.logWeights(particleIndex,interimDistIndex) = particleCollection.barLogWeights(particleIndex,interimDistIndex-1)+particleCollection.incremWeights(particleIndex,interimDistIndex);

      double particleAcceptanceRate =0;
      double propogatedParticles = mcmcMoveToy(particleCollection.samples(particleIndex, interimDistIndex-1),
                                               3,
                                               particleCollection.alphaVec(interimDistIndex),
                                               data,
                                               pow(stepSize,0.5),
                                               priorSigma, mu, likeSigma);
                                               // particleAcceptanceRate,
      particleCollection.accRate(interimDistIndex) = particleCollection.accRate(interimDistIndex)+particleAcceptanceRate;
      // std::cout<< "\n MCMC Kernel Completed";
      particleCollection.samples(arma::span(particleIndex,particleIndex),arma::span(interimDistIndex,interimDistIndex)) = propogatedParticles;
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

//'@export
// [[Rcpp::export]]
List smcSamplerToyCpp(int numParticles,int numOfDistrubtions, double dataSet,
                      int dimParameters, double priorSigma, double likeSigma,
                      double mu){

  SMCToyParticlesSystem  testingSMC = smcToy(numParticles, numOfDistrubtions, dataSet,
                                          dimParameters, priorSigma,
                                          likeSigma, mu);
  double normalizingConstant = arma::mean(exp(testingSMC.incremWeights.col(0)));
    for(int i = 1; i<numOfDistrubtions;i++){
      normalizingConstant = normalizingConstant * arma::sum(exp(testingSMC.incremWeights.col(i)+testingSMC.barLogWeights.col(i-1)));
    }
    // return(normalizingConstant);
  return(List::create(Named("particles") = testingSMC.samples, Named("normLogWeights") = testingSMC.logWeights, Named("UnNormlogWeights") = testingSMC.unNormLogWeights,
  Named("alphaVec") = testingSMC.alphaVec, Named("acceptanceRates") = testingSMC.accRate, Named("incremWeights") = testingSMC.incremWeights,
  Named("barLogWeights") = testingSMC.barLogWeights, Named("ESS") = testingSMC.ESS, Named("normalisingConstant") = normalizingConstant,
  Named("reSamplingIndex") = testingSMC.reSamplingIndices));
}

// int main(){
  // int arg =0 ;
  // double test = smcSamplerRcpp(200,500, data,timeFrames, 5);
  // cout<<test;


  // return(arg);
// }


// this takes the intiial sample generated using R -- since the sampler couldn't be found in C++
// [[Rcpp::export]]
void rTakeitToy(arma::vec test){
  // cout<< "\n rTake it :"<<rToySamples;
  rToySamples  = test;
}
