// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <string>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <cmath>
#include "truncated_normal.hpp"
#include "sourceVerOne.hpp"
#include <Rmath.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace std;
using namespace Rcpp;

int convMultiplier = 100001;
// float convMultiplier = 100001;
// arma::vec timeFrames;
// arma::mat convMatrix;
// bool success = convMatrix.load("convMatrix.Mat",arma::arma_ascii);
arma::vec convVector(100001*33);
bool successvec = convVector.load("./src/convVector2.vec",arma::raw_binary);

// arma::vec timeFrames = {0.0,27.5,60.0,70.0,80.0,100.0,130.0,160.0,190.0,220.0,250.0,280.0,355.0,
//               475.0,595.0,715.0,835.0,955.0,1075.0,1195.0,1315.0,1435.0,1555.0,1675.0,1885.0,2185.0,
//               2485.0,2785.0,3085.0,3385.0,3835.0,4435.0,5035.0};

// arma::vec data = {0.0008199750,0.0017871899 , 0.0018674304, -0.0006880768, -0.0004680085,  3.4993148899,  7.0149475981,  8.3965827114 , 9.0373437567 ,
//                   9.2884931257 ,9.3662051640  ,9.3268443765  ,9.4911005701,  9.1717982210,  8.8833974168,  8.5500954157,  8.3578397516,  8.1274822345,
//                   7.6957669549 , 7.4726012098  ,7.1896656682 , 6.8339554226,  6.5436026896  ,6.0752194184,  5.2872657085,  4.7246544526  ,4.1290334444,
//                   3.6299445843  ,3.1518110475,  2.6224344293,  1.9760028477,  1.5049115747};

double alpha = 1e-3, hpBeta = 1e-3;

// void matrixToVec(){
//
//   int index;
//   // std::cout<< "Matrix to Vector Function called" << std::endl;
//   // std::cout<< "n_rows:"<< convMatrix.n_rows << "\n"<< "n_cols:"<< convMatrix.n_cols;
//   // std::cout<< "vecotor size:"<< convVector.size();
//   for(int timeIndex = 0;timeIndex<33;timeIndex++){
//     for(int j =0;j<100001;j++){
//
//       index = (j*33)+timeIndex;
//       // std::cout<< index << "i:" << i << "j:" << j <<"\n";
//       convVector(index) = convMatrix(j,timeIndex);
//
//     }
//   }
//
// }

// This needs to be fixed
double convMatrixElement(float theta, int i){
  return convVector[(theta*33)+i];
}

double intpolFitValueCpp(int index, arma::vec phi, arma::vec theta){

  double fittedValue =0;
  // int index = *pI;
  int numberOfCompartments = theta.size();

  for(int d = 0; d<numberOfCompartments;d++){

    double thetaComponent = theta(d) * convMultiplier;
    int thetaLowerBound = floor(thetaComponent)+1;
    int thetaUpperBound = thetaLowerBound+1;
    double lowerMultiplier = phi(d) * (thetaComponent-thetaLowerBound);
    double upperMultiplie =  phi(d) * (thetaUpperBound-thetaComponent);

    // std::cout<< "thetaUpperBound: "<<thetaUpperBound <<"\n";
    // std::cout<< "*index "<<index <<"\n";
    // std::cout<< "convMatrixElement(thetaUpperBound,i): "<< convMatrixElement(thetaUpperBound,index) <<"\n";
    fittedValue = fittedValue + (lowerMultiplier*convMatrixElement(thetaUpperBound,index) + upperMultiplie*convMatrixElement(thetaLowerBound,index));

    // std::cout<< "fittedValue: "<< fittedValue <<"\n";
  }

  // std::cout<< "fittedValue: "<< fittedValue <<"\n";
  return(fittedValue);

}


double vaguePriorPhi(arma::vec phi){
  double priorDensity = 1;
  // have to use a different way here from R.
  for(int i = 0; i<phi.n_elem;i++){
    if( phi(i)<1e-5 || phi(i)>1 ){
      priorDensity = 0;
    }
  }
  return priorDensity;
}

double vaguePriorTheta(arma::vec theta){
  double priorDensity = 1;
  // have to use a different way here from R.
  for(int i = 0; i<theta.n_elem;i++){
    if( theta(i)<1e-4 || theta(i)>1 ){
      priorDensity = 0;
    }
  }
  return priorDensity;
}

double vaguePriorPhi(double phi){
  double priorDensity = 0;
  // have to use a different way here from R.
  if( phi>=1e-5 && phi<=1 ){
    priorDensity = 1;
  }
  return priorDensity;
}

double vaguePriorTheta(double theta){
  double priorDensity = 0;
  // have to use a different way here from R.
  if( theta>=1e-4 && theta<=1 ){
    priorDensity = 1;
  }
  return priorDensity;
}


double densityBioPriors(arma::vec phi, arma::vec theta){


  int dim  = phi.n_elem;
  for(int i = 0; i<dim;i++){
    if(phi(i)<0 || theta(i)<0){
      return(0);
    }
  }
  // if(vaguePriorTheta(theta) == 0 || vaguePriorPhi(phi) == 0){
  // return(0);
  // }
  double density;

  arma::vec phiPriors(dim,arma::fill::zeros);
  arma::vec thetaPriors(dim,arma::fill::zeros);

  double phiMean = 3e-3;
  double sigma = sqrt(1e-2);

  for(int i=0;i<dim;i++){
    phiPriors(i) = truncated_normal_ab_pdf(phi(i),phiMean,sqrt(1e-3),1e-5,1e-2);
    if(i==0){
      thetaPriors(i) = truncated_normal_ab_pdf(theta(i),(phi(i)/15),sigma,2e-4,1e-2);
      phiMean = 1e-3;
    }
    else if(i==1){
      thetaPriors(i) = truncated_normal_ab_pdf(theta(i),(phi(i)/4),sigma,theta(0),6e-2);
    }
    else if(i ==2){
      thetaPriors(i) = truncated_normal_ab_pdf(theta(i),(phi(i)),sigma,theta(i-1),6e-2);
    }
  }
  // std::cout<< "\n Theta Priors" << thetaPriors;
  // std::cout<< "\n Phi Priors" << phiPriors;

  density = arma::prod(thetaPriors) * arma::prod(phiPriors);
  return(density);
}

double petLogBioPriorDensity(arma::vec proposalState, int dim){
  int modelOrder = (dim-1)/2;
  arma::vec phi  = proposalState(arma::span(0,modelOrder-1));
  arma::vec theta = proposalState(arma::span(modelOrder,2*modelOrder-1));
  double b = proposalState(dim-1);
  // if(vaguePriorTheta(theta) == 0 || vaguePriorPhi(phi) == 0){
  // return(log(0));
  // }

  if(b <= -745.1332191 ){
    return(log(0));}

  for(int m = 0; m<modelOrder;m++){
    if(vaguePriorPhi(phi(m))==0 ||  vaguePriorTheta(theta(m)) ==0){
      // std::cout<<"\n vaguePriorPhi(phi(m)):"<<vaguePriorPhi(phi(m))<<" vaguePriorTheta(theta(m)):"<< vaguePriorTheta(theta(m));
      // std::cout<<"\n outsideVague: phi:"<<phi(m)<<" theta:"<< theta(m);
      return(log(0));
    }
  }
  return( (alpha*b)-(hpBeta*exp(b)) + log(densityBioPriors(phi,theta)));
}
// this only works with the look-up table optimisation, so no timeFrames arguement here.
double logLikelihood(arma::vec data, arma::vec timeFrames,arma::vec phi, arma::vec theta, double b){

  // we will need a summand as we are looking at the log transformed likelihood.
  double logLikelihood = 0;
  for(int i = 0; i < data.size();i++){

    double tissueConcentrationAtInterval = intpolFitValueCpp(i+1,phi,theta);
    // std::cout<< "tissueConcentrationAtInterval: "<< tissueConcentrationAtInterval <<"\n";
    tissueConcentrationAtInterval = fabs(tissueConcentrationAtInterval);

    logLikelihood = logLikelihood
      - 0.5*exp(b)*((timeFrames[i+1] - timeFrames[i])/tissueConcentrationAtInterval)*pow((data[i]-tissueConcentrationAtInterval),2)
      + 0.5*log((timeFrames[i+1] - timeFrames[i])/tissueConcentrationAtInterval);
      // std::cout<< "logLikelihood: "<< logLikelihood <<"\n";
  }


  logLikelihood = logLikelihood + data.size()*(b-log(2*3.141592653589793))*(0.5);

  return(logLikelihood);

}


double logPosterior(arma::vec data, arma::vec timeFrames, arma::vec phi, arma::vec theta, double b){

  double logPosteriorDensity  = 0;

  double logLikelihoodDensity = logLikelihood(data, timeFrames,phi,theta,b);

  logPosteriorDensity = logLikelihoodDensity + (alpha*b) -(hpBeta*exp(b));

  logPosteriorDensity  = logPosteriorDensity+log(densityBioPriors(phi,theta));

  return(logPosteriorDensity);

}

double logPosteriorWithAnnealing(arma::vec data, arma::vec timeFrames, arma::vec phi, arma::vec theta, double b, double alphaT){

  double logPosteriorDensity  = 0;

  double logLikelihoodDensity = logLikelihood(data, timeFrames,phi,theta,b);

  logPosteriorDensity = (alphaT*logLikelihoodDensity) + (alpha*b) -(hpBeta*exp(b));

  logPosteriorDensity  = logPosteriorDensity+log(densityBioPriors(phi,theta));

  return(logPosteriorDensity);

}
//
// int main(){
//
//
//   // arma::mat testSave;
//   // testSave.zeros(2,2);
//   // testSave.save("testing.MAT",arma::arma_ascii);
//
//   // matrixToVec();
//
//   std::cout<< "convMatrix:" <<convMatrix(1121,1) << "\n";
//
//
//
//   // int index  = (1*convMultiplier)+1;
//
//   // std::cout<< "convVectors1023:" <<convVector(10002*33+10);
//   // convVector.save("convVector.vec",arma::arma_ascii);
//
//
//   // float b  = 2.538946;
//
//   // arma::vec phi = {0.0032, 0.000655};
//   // arma::vec theta  ={0.054654, 0.2132} ;
//
//   // arma::vec phi = {0.1};
//   // arma::vec theta  ={0.1} ;
//
//
//   // int ten = 10;
//   // std::cout<< "LogLikelihood:"<< logLikelihood(data,timeFrames,phi,theta,b) << "\n";
//   // std::cout<< "\n fitted Value:" << intpolFitValue(&ten,phi,theta);
//   // std::cout<< "\n This is 1e-2" << 1e-2;
//   // std::cout<< "LogLikelihood:"<< logLikelihood(data,timeFrames,phi,theta,-2)<< "\n";
//   // std::cout<< "LogLikelihood:"<< logLikelihood(data,timeFrames,phi,theta,-3)<< "\n";
//   // std::cout<< "LogLikelihood:"<< logLikelihood(data,timeFrames,phi,theta,-1)<< "\n";
//   // std::cout<< "LogLikelihood:"<< logLikelihood(data,timeFrames,phi,theta,-0.3) << "\n";
//
//
//   // densityBPriors testing
//   // std::cout<< "\n Denisty of Bio Priors for phi: "<<phi<< "and theta: "<<theta<<  " is :"<<densityBioPriors(phi,theta);
//   arma:: vec phi  = {0.01};
//   arma::vec theta  = {0.01};
//   // std::cout<< "\n Denisty of Bio Priors for phi: "<<phi<< "and theta: "<<theta<<  " is :"<<densityBioPriors(phi,theta);
//
//
//   std::cout<< "\n Log Posterior Density: "<<phi<< "and theta: "<<theta<<  " is :"<<logPosterior(data,timeFrames,phi,theta,log(0.64));
//
//  //  phi  = {0.01};
//  //  theta  = {0.01};
//  //  std::cout<< "\n Denisty of Bio Priors for phi: "<<phi<< "and theta: "<<theta<<  " is :"<<densityBioPriors(phi,theta);
//  //  phi  = {0.01,0.01};
//  //  theta  = {0.01,0.01};
//  //  std::cout<< "\n Denisty of Bio Priors for phi: "<<phi<< "and theta: "<<theta<<  " is :"<<densityBioPriors(phi,theta);
//  //  phi  = {0.01,0.01,0.01};
//  //  theta  = {0.01,0.01,0.1};
//  //  std::cout<< "\n Denisty of Bio Priors for phi: "<<phi<< "and theta: "<<theta<<  " is :"<<densityBioPriors(phi,theta);
//  //  phi  = {0.01,0.05,0.02};
//  //  theta  = {0.01,0.02,0.05};
//  //  std::cout<< "\n Denisty of Bio Priors for phi: "<<phi<< "and theta: "<<theta<<  " is :"<<densityBioPriors(phi,theta);
//  std::cout<< "\n HELLO";
//   return 0;
// }


// arma::vec csample(int numOfSamples, arma::vec prob){
//   std::vector<double> tempProb = arma::conv_to<std::vector<double>>::from(prob);
//   NumericVector probNumVec = wrap(tempProb);
//   IntegerVector indices = seq(0,numOfSamples-1);
//   IntegerVector ret = sample(indices,numOfSamples,true,probNumVec);
//
//
//   std::vector<double> tempRet = as<std::vector<double> >(ret);
//   arma::vec realRet = arma::conv_to<arma::vec>::from(tempRet);
//
//   return realRet;
// }

IntegerVector csampleNoConv(IntegerVector indices, NumericVector prob){
  IntegerVector ret = sample(indices,10,true,prob);
  return ret;
}

arma::uvec csample(int numOfSamples, arma::vec prob){
  std::vector<double> tempProb = arma::conv_to<std::vector<double>>::from(prob);
  NumericVector probNumVec = wrap(tempProb);
  IntegerVector indices = seq(0,numOfSamples-1);
  IntegerVector ret = sample(indices,numOfSamples,true,probNumVec);


  std::vector<double> tempRet = as<std::vector<double> >(ret);
  arma::uvec realRet = arma::conv_to<arma::uvec>::from(tempRet);

  return realRet;
}
