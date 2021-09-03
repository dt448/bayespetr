// [[Rcpp::depends(RcppArmadillo)]]
//============================================================================
// vSMC/example/pet/src/pet_smc.cpp
//----------------------------------------------------------------------------
//                         vSMC: Scalable Monte Carlo
//----------------------------------------------------------------------------
// Copyright (c) 2013-2015, Yan Zhou
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//   Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
//
//   Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//============================================================================

#include "pet_seq.hpp"
#include "smc.hpp"
// #include <stdlib.h>     /* srand, rand */
// #include <time.h>       /* time */
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

// std::string FILEPATH = "./";

std::vector<double> pet_data_from_r;
std::vector<double> time_interv_from_r;
std::vector<double> prior_ranges_from_r;
double decay_low_rate_from_r;

// std::vector<double> Time{
// 27.5, 32.5, 10.0, 10.0, 20.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 75.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0,
// 120.0, 120.0, 120.0, 120.0, 120.0, 210.0, 300.0, 300.0, 300.0, 300.0, 300.0, 450.0, 600.0, 600.0};

// std::vector<double> Time{
//   40, 45, 10.0, 10.0, 20.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 75.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0,
//   120.0, 120.0, 120.0, 120.0, 120.0, 210.0, 300.0, 300.0, 300.0, 300.0, 300.0, 450.0, 600.0, 600.0};


// std::vector<double> Time{ 80,  90, 100,110, 140, 170, 200, 230, 260, 290, 320,
//                           440, 560, 680, 800, 920, 1040, 1160, 1280, 1400, 1520,
//                           1640, 1760, 2060, 2360, 2660, 2960, 3260, 3560, 4160,
//                           4760, 5360};


//h01031_timeFrames:
// std::vector<double> Time{ 45, 55,65,75,105,135,165,195,225,255,285,405,525,645,765,
//                           885,1005,1125,1245,1365,1485,1605,1725,2025,2325,
//                           2625,2925,3225,3525,4125,4725,5325};

std::vector<double> SD{ 5e-4, 1e-3, 1e-3,
                        5e-5, 1e-2, 1e-2,
                        0.3,
                        0.3};

// original given prior
// std::vector<double> Prior{1e-3, 1e-5, 1e-5,
//                           0.01, 0.1, 0.1,
//                           1e-4, 1e-5, 1e-5,
//                           1e-3, 1e-1, 1e-1,
//                           1e-1,
//                           10,
//                           0,
//                           5e-1};
// std::vector<double> Prior{1e-5, 1e-5, 1e-5,
//                           1, 1, 1,
//                           1e-4, 1e-4, 1e-4,
//                           1, 1, 1,
//                           1e-1,
//                           10,
//                           0,
//                           5e-1};

// std::vector<double> Prior{1e-5, 1e-5, 1e-5,
//                           0.01, 0.01, 0.01,
//                           1e-4, 1e-4, 1e-4,
//                           0.01, 0.06, 0.06,
//                           1e-1,
//                           10,
//                           0,
//                           5e-1};


void print(std::vector <double> const &a) {
  // std::cout << "The vector elements are : ";

  for(int i=0; i < a.size(); i++)
    std::cout << a.at(i) << ' ';
}

void print(std::vector <long unsigned int> &a) {
  // std::cout << "The vector elements are : ";

  for(int i=0; i < a.size(); i++)
    std::cout << a.at(i) << ' ';
}

//'@export
// [[Rcpp::export]]
void dataFromR(NumericVector r_pet_data, List config_model){
  pet_data_from_r = as<std::vector<double>>(r_pet_data);
  time_interv_from_r = as<std::vector<double>>(config_model["time_intervals"]);
  prior_ranges_from_r = as<std::vector<double>>(config_model["prior_ranges"]);
  decay_low_rate_from_r = config_model["decay_lowest_rate"];
  // print(DataFromR);
}

std::vector<double> readConvFile(std::string PETCONV_FILEPATH){

  std::size_t ConvNum = 100002;
  std::size_t DataNum = 32;
  std::vector<double> Conv(32*100002);
  std::ifstream data_file;
  data_file.open(PETCONV_FILEPATH);
  // data_file.open(FILEPATH+"pet_conv_h02913.data");
  // data_file.open(FILEPATH+"new_pet_conv_2.data");
  // data_file.open(FILEPATH+"pet_conv_h21903_timeCorrected.data");
  // data_file.open(FILEPATH+"pet_conv_time_corrected_end_tfs.data");
  // data_file.open(FILEPATH+"h01031_pet_conv.data");

  for (std::size_t r = 0; r != ConvNum; ++r)
    for (std::size_t c = 0; c != DataNum; ++c)
      data_file >> Conv[c + r * DataNum];
  data_file.close();
  data_file.clear();
  return Conv;
}


std::vector<double> ConvPass;

//'@export
// [[Rcpp::export]]
void cpp_load_conv_table(std::string pet_conv_dir){
  ConvPass = readConvFile(pet_conv_dir);
}

//'@export
// [[Rcpp::export]]
List vSMCRCppSampler(std::size_t numSamples, double prior5_iter_num, std::size_t modelNumber, List R_config, int randSeed){

  std::string MCMC_kernel_type = R_config["mcmc_type"];
  std::string annealing_scheme = R_config["annealing_scheme"];
  double cess_threshold;
  int mcmc_iters =  R_config["mcmc_iters"];
  // vsmc takes in tuning parameters via the Config(ProgrammMap Option)
  // to keep things consistent (i.e. prevent having to chnage vSMC itself)
  // the following will parse everything over.
  if( annealing_scheme == "CESS"){
    cess_threshold = R_config["cess_threshold"];
    CESSDrop.clear();
    CESSDrop.push_back(1.0-cess_threshold);
  }else if(annealing_scheme == "User"){
    alpha_user_input = as<std::vector<double>>(R_config["user_annealing_scheme"]);
  }else{
    // slightly longer but will sanatize arguements.
    // std::size_t prior5_iter_num = R_config["prior5_iter_num"];
    Prior5IterNum.clear();
    Prior5IterNum.push_back(prior5_iter_num);
  }

  // TODO : may not need this
#include "options_main.hpp"
  vsmc::Seed::instance().set(randSeed);
  // std::cout<< "importing options_main.hpp completed"<<endl;
  // Prior5IterNum.clear();
  // Prior5IterNum.push_back(numDistributions);
  // std::cout<<"Prior after pushback"<<std::endl;
  // print(Prior5IterNum);
  ParticleNum = numSamples;
#include "options_smc.hpp"
  // std::cout<< "importing options_smc.hpp completed"<<endl;
  // #include "pet_options.hpp"
  // DataFile = "/home/den/Desktop/dataFolder/pet_simulate_data.config";
  // DataFile = FILEPATH+"pet_simulate_data.config";
  SM = 1;
  CM = 2;
  // std::cout<< "importing pet_options.hpp completed"<<endl;
#include "options_process.hpp"
  // std::cout<< "importing options_process.hpp completed"<<endl;
  // std::cout<< "importing pet_data.hpp completed \n";
  // std::cout << "pet_data.hpp DataFile" <<DataFile<< std::endl;
#include "pet_data_Rcpp.hpp"
  // passing model information into the pet_info struct
  std::vector<double> Data = pet_data_from_r;
  data_info  info_d = {DataNum, &Data[0]};

  std::vector<double> Time = time_interv_from_r;
  time_info  info_t = {DataNum, &Time[0]};;

  std::vector<double> Conv = ConvPass;
  conv_info  info_c = {ConvNum, ConvMul, &Conv[0]};;

  std::vector<double> Prior = prior_ranges_from_r;
  prior_info info_p = {ModelNum, &Prior[0]};;
  Decay = decay_low_rate_from_r;

  PETModel ModelType;

  if( R_config["t_distribution"])
  {
    ModelType = StudentT;
  }else{
    ModelType = Normal;
  }
  model_info info_m = {Decay, ModelType};;
  pet_info info = {
    &info_d, true,
    &info_t, true,
    &info_c, true,
    &info_p, true,
    &info_s, true,
    &info_m, true
  };

  double normalizingConstant;
  double volume_of_distribution;

  arma::vec phiOutput(modelNumber);
  arma::vec thetaOutput(modelNumber);

  if(MCMC_kernel_type.compare("AMCMC") == 0){
    ProposalScale = 2;
  }else{
    ProposalScale = 1;
  }

  vsmc::Sampler<pet_state> sampler(ParticleNum, vsmc::Stratified, Threshold);
  sampler
    .init(pet_init())
    .mcmc(pet_move_phi(), true)
    .mcmc(pet_move_theta(), true)
    .mcmc(pet_move_lambda(), true)
    .mcmc(pet_move_nu(), true);
  if(mcmc_iters>1){
    for(int i = 1; i <mcmc_iters;++i){
      sampler
      .mcmc(pet_move_phi(), true)
      .mcmc(pet_move_theta(), true)
      .mcmc(pet_move_lambda(), true)
      .mcmc(pet_move_nu(), true);
    }
  }


  sampler
    .monitor("vd", 2, pet_vd())
    .monitor("phi", modelNumber, pet_phi())
    .monitor("theta", modelNumber, pet_theta())
    .path_sampling(smc_path<pet_state>());

  if (ProposalScale == 2) {
    sampler.monitor("pet_moments", 2 * (2 + 2 * modelNumber),
                    pet_moments());
  }
  //
  // std::cout<<"\n Initializing done !";
  // std::cout<<"\n InitCompNum"<<InitCompNum;
  // std::cout<<"\n InitCompNumBARNEW "<<InitCompNumBAR;
  sampler.initialize(&info);
  info.read_time  = false;
  info.read_conv  = false;
  info.read_prior = false;
  info.read_sd    = false;
  info.read_model = false;

  // std::cout<<"\n Initializing initiing done !";
  //
  //     //////////////////////////////////////////////////////////////////////
  //
  std::string zconst_file_name("~/smc." + Suffix);
  std::ofstream zconst_file;
  zconst_file.open(zconst_file_name.c_str());
  zconst_file << "Schedule Config ";
  print_zconst_header(zconst_file, modelNumber);
  // print_zconst_header(zconst_file, CM);
  zconst_file << std::endl;

  if (ProposalScale == 2) {
    typedef pet_proposal_adaptive sd;
    for (std::size_t i = DataStart; i != DataStop; ++i) {
      info_d.data_value = &Data[i * DataNum];
      sampler.initialize(&info);
      smc_do<pet_state, sd>(annealing_scheme, sampler, zconst_file, modelNumber);
    }
  } else{
    typedef pet_proposal sd;
    for (std::size_t i = DataStart; i != DataStop; ++i) {
      info_d.data_value = &Data[i * DataNum];
      sampler.initialize(&info);

      //             // sampler.iterate(100);
      //             // vsmc::Particle<pet_state> &particle = sampler.particle();
      //             // double average = sampler.log_zconst();
      //             // std::cout<<average;
      //             // std::size_t modelNumber = 2;

      smc_do<pet_state, sd>(annealing_scheme, sampler, zconst_file, modelNumber);

      // std::cout<<"normalizingConstant:"<<normalizingConstant;
      //             // smc_do<pet_state>(sampler, zconst_file,"test",3,500);
      //
      //
      //
      //             // double  zconst = sampler.particle().value().zconst();
      //             // double  pathLogzconst = sampler.path().log_zconst();
      //             // std::size_t iterSize = sampler.iter_size();
      //             // std::cout<<"zconst:"<<zconst;
      //             // std::cout<<"log_likelihood_const:"<<sampler.particle().value().log_likelihood_const();
      //             // std::cout<<"iterSize:"<<iterSize;
      //             // double *lw = new double[particle.size()];
      //             // particle.weight_set().read_log_weight(lw);
      //             // std::cout<<"logweight[1]"<<lw[0];
      //
      //             // std::vector<double> w(particle.size());
      //             // particle.weight_set().read_weight(w.begin());
      //             // std::cout<<"weight[1]"<<w[0];
      //
      //
    }
  }

  normalizingConstant = sampler.particle().value().zconst()
    + sampler.particle().value().log_likelihood_const();
  volume_of_distribution = sampler.monitor("vd").record(0);

  for(int i = 0; i<modelNumber; i++){
    phiOutput[i] = sampler.monitor("phi").record(i);
    thetaOutput[i] = sampler.monitor("theta").record(i);
    // std::cout<<"phi0:"<<phiOutput[i]<<"theta0:"<<thetaOutput[i];
  }

  zconst_file.close();
  zconst_file.clear();
  return(List::create(Named("normalizingConstant")=normalizingConstant,
                      Named("volumeOfDistribution") = volume_of_distribution,
                      Named("phi") = phiOutput,
                      Named("theta") = thetaOutput,
                      Named("annealing_scheme") = alpha_vec
  ));
}
