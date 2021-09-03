#include <iostream>
#include <RcppArmadillo.h>
#include "sourceVerOne.hpp"

arma::vec convVector2(100002*32);
// bool successvec = convVector.load("/home/u1760857/bayespetr/convVector2.vec",arma::raw_binary);
bool successvec2 = convVector2.load("/home/den/Desktop/dataFolder/convVector_h01031.vec",arma::raw_binary);

double convMatrixElementFixedRcpp(float theta, int i){
  return convVector2[(theta*32)+i];
}

// [[Rcpp::export]]
double intpolFitValue(Rcpp::IntegerVector timeIndex, Rcpp::NumericVector phi, Rcpp::NumericVector theta){

//************************************WARNING:Any changes made here should also
//**************be made in the sourceVerOne.cpp file****************************
  double fittedValue =0;
  int index = timeIndex[0] -1;
  int numberOfCompartments = theta.length();

  for(int d = 0; d<numberOfCompartments;d++){

    double thetaComponent = theta(d) * 100000;
    int thetaLowerBound = floor(thetaComponent);
    int thetaUpperBound = thetaLowerBound+1;
    double lowerMultiplier = phi(d) * (thetaComponent-thetaLowerBound);
    double upperMultiplie =  phi(d) * (thetaUpperBound-thetaComponent);


    // cout.precision(dbl::max_digits10);
    // std::cout<< "thetaUpperBound: "<<thetaUpperBound <<"\n";
    // std::cout<< "upperMultiplier: "<<upperMultiplie <<"\n";
    // std::cout<< "lowerMultiplier: "<<lowerMultiplier <<"\n";
    // std::cout<< "*index "<<index <<"\n";
    // std::cout<< "convMatrixElement(thetaUpperBound,i): "<< convMatrixElement(thetaUpperBound,index) <<"\n";
    // fittedValue = fittedValue + (lowerMultiplier*convMatrixElement(thetaUpperBound,index) + upperMultiplie*convMatrixElement(thetaLowerBound,index));
    fittedValue = fittedValue + (lowerMultiplier*convMatrixElementFixedRcpp(thetaUpperBound,index) + upperMultiplie*convMatrixElementFixedRcpp(thetaLowerBound,index));

    // std::cout<< "fittedValue: "<< fittedValue <<"\n";
  }

  // std::cout<< "fittedValue: "<< fittedValue <<"\n";
  return(fittedValue);

}


// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}
