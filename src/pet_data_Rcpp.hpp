//============================================================================
// vSMC/example/pet/include/pet_data.hpp
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

// std::cout << "pet_data_Rcpp.hpp called" << std::endl;


// Uses the external files just for decay, model num, data stop, datastart etc, clean this up and should work without needing
// external files



double Decay = 0;
std::size_t DataSetNum = 1; //dont need
std::size_t DataStart = 0; // dont need
std::size_t DataStop = 1; //dont need
std::size_t DataNum = 32; //base this on timeframes
std::size_t ModelNum = 3; // done need
std::size_t ConvNum = 100002;
double ConvMul = 100000;
// std::cout << "pet_data.hpp before strings" << std::endl;
/*std::string TimeFile = 0;
std::string ConvFile = 0;
std::string PriorFile = 0;
std::string SDFile = 0;*/
std::string TimeFile;
std::string ConvFile;
std::string PriorFile;
std::string SDFile;
bool UseStudentT = false;
// std::cout << "pet_data.hpp after strings" << std::endl;

// std::ifstream config_file;
// std::cout << "pet_data.hpp before opening DataFile" << std::endl;
// config_file.open(DataFile.c_str());
// std::cout << "pet_data.hpp After Datafile" << std::endl;
// config_file >> Decay >> DataSetNum >> DataStart >> DataStop
    // >> ModelNum >> DataNum >> ConvNum >> ConvMul
    // >> DataFile >> TimeFile >> ConvFile >> PriorFile >> SDFile >> UseStudentT;
// config_file.close();
// config_file.clear();

// std::cout << "pet_data.hpp DataFile" <<DataFile<< std::endl;
// std::cout<< "PriorFile:"<<PriorFile;


// std::vector<double> Time(DataNum);
// std::vector<double> Time{
//     27.5, 32.5, 10, 10, 20, 30, 30, 30, 30, 30, 30, 75, 120, 120, 120, 120, 120, 120,
//     120, 120, 120, 120, 120, 210, 300, 300, 300, 300, 300, 450, 600, 600};
// // std::vector<double> Conv(DataNum * ConvNum);
//
// std::vector<double> Prior{1e-3, 1e-5, 1e-5,
//                           1e-2, 1e-1, 1e-1,
//                           1e-4, 1e-5, 1e-5,
//                           1e-3, 1e-1, 1e-1,
//                           1e-1,
//                           10,
//                           0,
//                           5e-1};
//
//
//
// // std::vector<double> SD(ModelNum * 2 + 2);
// std::vector<double> SD{ 5e-4, 1e-3, 1e-3,
//                         5e-5, 1e-2, 1e-2,
//                         0.3,
//                         0.3};
// if (UseStudentT)
    // ModelType = StudentT;

// std::vector<double> Data(DataNum * DataSetNum);
//
// data_file.open(DataFile.c_str());
// for (std::size_t r = 0; r != DataSetNum; ++r)
//     for (std::size_t c = 0; c != DataNum; ++c)
//         data_file >> Data[c + r * DataNum];
// data_file.close();
// data_file.clear();

// data_file.open(TimeFile.c_str());

// std::ifstream data_file;
// for (std::size_t i = 0; i != DataNum; ++i)
    // data_file >> Time[i];
// data_file.close();
// data_file.clear();

// data_file.open(ConvFile.c_str());
// for (std::size_t r = 0; r != ConvNum; ++r)
    // for (std::size_t c = 0; c != DataNum; ++c)
        // data_file >> Conv[c + r * DataNum];
// data_file.close();
// data_file.clear();

// std::size_t prior_offset = 0;
// data_file.open(PriorFile.c_str());
// for (std::size_t i = 0; i != ModelNum; ++i)
    // data_file >> Prior[prior_offset++]; pet_ignore(data_file); // phi_lb0
// for (std::size_t i = 0; i != ModelNum; ++i)
    // data_file >> Prior[prior_offset++]; pet_ignore(data_file); // phi_ub0
// for (std::size_t i = 0; i != ModelNum; ++i)
    // data_file >> Prior[prior_offset++]; pet_ignore(data_file); // theta_lb0
// for (std::size_t i = 0; i != ModelNum; ++i)
    // data_file >> Prior[prior_offset++]; pet_ignore(data_file); // theta_ub0
// data_file >> Prior[prior_offset++]; pet_ignore(data_file); // lambda_a0
// data_file >> Prior[prior_offset++]; pet_ignore(data_file); // lambda_b0
// data_file >> Prior[prior_offset++]; pet_ignore(data_file); // nu_a0
// data_file >> Prior[prior_offset++]; pet_ignore(data_file); // nu_b0
// data_file.close();
// data_file.clear();

// std::size_t sd_offset = 0;
// data_file.open(SDFile.c_str());
// for (std::size_t i = 0; i != ModelNum; ++i)
    // data_file >> SD[sd_offset++]; pet_ignore(data_file); // phi_sd
// for (std::size_t i = 0; i != ModelNum; ++i)
    // data_file >> SD[sd_offset++]; pet_ignore(data_file); // theta_sd
// data_file >> SD[sd_offset++]; pet_ignore(data_file); // lamda_sd;
// data_file >> SD[sd_offset++]; pet_ignore(data_file); // nu_sd;
// data_file.close();
// data_file.clear();


sd_info    info_s = {ModelNum, &SD[0]};;
// pet_info info = {
//     &info_d, true,
//     &info_t, true,
//     &info_c, true,
//     &info_p, true,
//     &info_s, true,
//     &info_m, true
// };

