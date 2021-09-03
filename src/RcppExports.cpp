// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// smcSamplerCpp
List smcSamplerCpp(int numParticles, int numOfDistrubtions, arma::vec dataSet, arma::vec timeFrames, int dimParameters);
RcppExport SEXP _bayespetr_smcSamplerCpp(SEXP numParticlesSEXP, SEXP numOfDistrubtionsSEXP, SEXP dataSetSEXP, SEXP timeFramesSEXP, SEXP dimParametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numParticles(numParticlesSEXP);
    Rcpp::traits::input_parameter< int >::type numOfDistrubtions(numOfDistrubtionsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dataSet(dataSetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type timeFrames(timeFramesSEXP);
    Rcpp::traits::input_parameter< int >::type dimParameters(dimParametersSEXP);
    rcpp_result_gen = Rcpp::wrap(smcSamplerCpp(numParticles, numOfDistrubtions, dataSet, timeFrames, dimParameters));
    return rcpp_result_gen;
END_RCPP
}
// rTakeit
void rTakeit(arma::mat test);
RcppExport SEXP _bayespetr_rTakeit(SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type test(testSEXP);
    rTakeit(test);
    return R_NilValue;
END_RCPP
}
// mcmcMoveToy
double mcmcMoveToy(double initialValue, int steps, double alphaT, double data, double stepSize, double priorSigma, double mu, double likeSigma);
RcppExport SEXP _bayespetr_mcmcMoveToy(SEXP initialValueSEXP, SEXP stepsSEXP, SEXP alphaTSEXP, SEXP dataSEXP, SEXP stepSizeSEXP, SEXP priorSigmaSEXP, SEXP muSEXP, SEXP likeSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type initialValue(initialValueSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type alphaT(alphaTSEXP);
    Rcpp::traits::input_parameter< double >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type stepSize(stepSizeSEXP);
    Rcpp::traits::input_parameter< double >::type priorSigma(priorSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type likeSigma(likeSigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmcMoveToy(initialValue, steps, alphaT, data, stepSize, priorSigma, mu, likeSigma));
    return rcpp_result_gen;
END_RCPP
}
// smcSamplerToyCpp
List smcSamplerToyCpp(int numParticles, int numOfDistrubtions, double dataSet, int dimParameters, double priorSigma, double likeSigma, double mu);
RcppExport SEXP _bayespetr_smcSamplerToyCpp(SEXP numParticlesSEXP, SEXP numOfDistrubtionsSEXP, SEXP dataSetSEXP, SEXP dimParametersSEXP, SEXP priorSigmaSEXP, SEXP likeSigmaSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numParticles(numParticlesSEXP);
    Rcpp::traits::input_parameter< int >::type numOfDistrubtions(numOfDistrubtionsSEXP);
    Rcpp::traits::input_parameter< double >::type dataSet(dataSetSEXP);
    Rcpp::traits::input_parameter< int >::type dimParameters(dimParametersSEXP);
    Rcpp::traits::input_parameter< double >::type priorSigma(priorSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type likeSigma(likeSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(smcSamplerToyCpp(numParticles, numOfDistrubtions, dataSet, dimParameters, priorSigma, likeSigma, mu));
    return rcpp_result_gen;
END_RCPP
}
// rTakeitToy
void rTakeitToy(arma::vec test);
RcppExport SEXP _bayespetr_rTakeitToy(SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type test(testSEXP);
    rTakeitToy(test);
    return R_NilValue;
END_RCPP
}
// intpolFitValue
double intpolFitValue(Rcpp::IntegerVector timeIndex, Rcpp::NumericVector phi, Rcpp::NumericVector theta);
RcppExport SEXP _bayespetr_intpolFitValue(SEXP timeIndexSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type timeIndex(timeIndexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(intpolFitValue(timeIndex, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// set_seed
void set_seed(unsigned int seed);
RcppExport SEXP _bayespetr_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// dataFromR
void dataFromR(NumericVector r_pet_data, List config_model);
RcppExport SEXP _bayespetr_dataFromR(SEXP r_pet_dataSEXP, SEXP config_modelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r_pet_data(r_pet_dataSEXP);
    Rcpp::traits::input_parameter< List >::type config_model(config_modelSEXP);
    dataFromR(r_pet_data, config_model);
    return R_NilValue;
END_RCPP
}
// cpp_load_conv_table
void cpp_load_conv_table(std::string pet_conv_dir);
RcppExport SEXP _bayespetr_cpp_load_conv_table(SEXP pet_conv_dirSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type pet_conv_dir(pet_conv_dirSEXP);
    cpp_load_conv_table(pet_conv_dir);
    return R_NilValue;
END_RCPP
}
// vSMCRCppSampler
List vSMCRCppSampler(std::size_t numSamples, double prior5_iter_num, std::size_t modelNumber, List R_config, int randSeed);
RcppExport SEXP _bayespetr_vSMCRCppSampler(SEXP numSamplesSEXP, SEXP prior5_iter_numSEXP, SEXP modelNumberSEXP, SEXP R_configSEXP, SEXP randSeedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::size_t >::type numSamples(numSamplesSEXP);
    Rcpp::traits::input_parameter< double >::type prior5_iter_num(prior5_iter_numSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type modelNumber(modelNumberSEXP);
    Rcpp::traits::input_parameter< List >::type R_config(R_configSEXP);
    Rcpp::traits::input_parameter< int >::type randSeed(randSeedSEXP);
    rcpp_result_gen = Rcpp::wrap(vSMCRCppSampler(numSamples, prior5_iter_num, modelNumber, R_config, randSeed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayespetr_smcSamplerCpp", (DL_FUNC) &_bayespetr_smcSamplerCpp, 5},
    {"_bayespetr_rTakeit", (DL_FUNC) &_bayespetr_rTakeit, 1},
    {"_bayespetr_mcmcMoveToy", (DL_FUNC) &_bayespetr_mcmcMoveToy, 8},
    {"_bayespetr_smcSamplerToyCpp", (DL_FUNC) &_bayespetr_smcSamplerToyCpp, 7},
    {"_bayespetr_rTakeitToy", (DL_FUNC) &_bayespetr_rTakeitToy, 1},
    {"_bayespetr_intpolFitValue", (DL_FUNC) &_bayespetr_intpolFitValue, 3},
    {"_bayespetr_set_seed", (DL_FUNC) &_bayespetr_set_seed, 1},
    {"_bayespetr_dataFromR", (DL_FUNC) &_bayespetr_dataFromR, 2},
    {"_bayespetr_cpp_load_conv_table", (DL_FUNC) &_bayespetr_cpp_load_conv_table, 1},
    {"_bayespetr_vSMCRCppSampler", (DL_FUNC) &_bayespetr_vSMCRCppSampler, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayespetr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}