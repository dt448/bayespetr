// [[Rcpp::depends(RcppArmadillo)]]
#ifndef FILE_FOO_SEEN
#define FILE_FOO_SEEN
#include <iostream>
#include <string>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <cmath>
#include "truncated_normal.hpp"
#include <Rmath.h>

// arma::vec timeFrames = {0.0,27.5,60.0,70.0,80.0,100.0,130.0,160.0,190.0,220.0,250.0,280.0,355.0,
//               475.0,595.0,715.0,835.0,955.0,1075.0,1195.0,1315.0,1435.0,1555.0,1675.0,1885.0,2185.0,
//               2485.0,2785.0,3085.0,3385.0,3835.0,4435.0,5035.0};
//
// arma::vec data = {0.0008199750,0.0017871899 , 0.0018674304, -0.0006880768, -0.0004680085,  3.4993148899,  7.0149475981,  8.3965827114 , 9.0373437567 ,
//                   9.2884931257 ,9.3662051640  ,9.3268443765  ,9.4911005701,  9.1717982210,  8.8833974168,  8.5500954157,  8.3578397516,  8.1274822345,
//                   7.6957669549 , 7.4726012098  ,7.1896656682 , 6.8339554226,  6.5436026896  ,6.0752194184,  5.2872657085,  4.7246544526  ,4.1290334444,
//                   3.6299445843  ,3.1518110475,  2.6224344293,  1.9760028477,  1.5049115747};

extern int convMultiplier;
extern double alpha, beta;
extern arma::vec convVector;
void matrixToVec();

double convMatrixElement(float theta, int i);

double intpolFitValue(int index, arma::vec phi, arma::vec theta);

double vaguePriorPhi(arma::vec phi);

double vaugePriorTheta(arma::vec theta);

double vaguePriorPhi(arma::vec phi);

double vaguePriorTheta(arma::vec theta);

double vaguePriorPhi(double phi);

double vaguePriorTheta(double theta);

double densityBioPriors(arma::vec phi, arma::vec theta);

double petLogBioPriorDensity(arma::vec proposalState, int dim);

// this only works with the look-up table optimisation, so no timeFrames arguement here.
double logLikelihood(arma::vec data, arma::vec timeFrames,arma::vec phi, arma::vec theta, double b);

double logPosterior(arma::vec data, arma::vec timeFrames, arma::vec phi, arma::vec theta, double b);

double logPosteriorWithAnnealing(arma::vec data, arma::vec timeFrames, arma::vec phi, arma::vec theta, double b, double alphaT);
arma::uvec csample(int numOfSamples, arma::vec prob);

#endif
