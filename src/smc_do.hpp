//============================================================================
// vSMC/example/include/smc_do.hpp
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

#ifndef VSMC_EXAMPLE_SMC_DO_HPP
#define VSMC_EXAMPLE_SMC_DO_HPP

#include "pet_smc_seq.hpp"
#include <RcppArmadillo.h>

using namespace Rcpp;

inline void print_zconst_header (std::ostream &zconst_file,
        std::size_t model_num)
{
    if (model_num == 0)
        return;

    zconst_file
        << "Iteration." << model_num
        << " SMC." << model_num << " Path." << model_num;
#ifdef VSMC_PET_HPP
    zconst_file << " VDMean." << model_num << " VDVar." << model_num;
#endif
    zconst_file << ' ';
}

template <typename T>
inline void smc_do (vsmc::Sampler<T> &sampler, std::ostream &zconst_file,
        const std::string schedule_name, std::size_t model_num,
        std::size_t iter_num)
{
    if (model_num == 0)
        return;

    sampler.particle().value().comp_num(model_num);
    sampler.initialize();
    // vsmc::StopWatch watch;
    // watch.start();
    alpha_vec.clear();
    if (iter_num == 0) {
        while (sampler.particle().value().state(0, 0).alpha() < 1){
            alpha_vec.push_back(sampler.particle().value().state(0, 0).alpha());
            sampler.iterate();
            Rcpp::checkUserInterrupt();
        }
    } else {
        sampler.iterate(iter_num);
    }
    // watch.stop();
    // std::fprintf(stderr, "time.model.order.%u %f\n",
            // static_cast<unsigned>(model_num), watch.seconds());
    // std::fflush(stderr);

    zconst_file << sampler.iter_size() - 1 << ' ';
    zconst_file << sampler.particle().value().zconst()
        + sampler.particle().value().log_likelihood_const() << ' ';

    zconst_file << sampler.path_sampling()
        + sampler.particle().value().log_likelihood_const() << ' ';
#ifdef VSMC_PET_HPP
    zconst_file << sampler.monitor("vd").record(0) << ' ';
    zconst_file << sampler.monitor("vd").record(1) -
        sampler.monitor("vd").record(0) * sampler.monitor("vd").record(0);
#endif
    zconst_file << ' ';

    // std::cout<<"zconst in side smc_do:"<<sampler.particle().value().zconst();
    // std::cout<<"log_likelihood_const:"<<sampler.particle().value().log_likelihood_const();

    // If you want the samples.....

    // std::ofstream sampler_file;
    // std::string fn;
    // std::stringstream ss;
    // ss << model_num;
    // std::string model_name(ss.str());

    // fn = std::string("smc.sampler.");
    // fn += Suffix + "." + schedule_name +  "." + model_name;
    // sampler_file.open(fn.c_str());
    // sampler_file << sampler;
    // sampler_file.close();
    // sampler_file.clear();
}

// modified from original to include model_num arguement
template <typename MoveType, typename T>
inline void smc_do (vsmc::Sampler<T> &sampler,
        std::ostream &zconst_file, const std::string &schedule_name,
        typename MoveType::alpha_type::value_type alpha_config,
        std::size_t iter_num, std::size_t model_num)
{
    sampler.move(MoveType(alpha_config, &sampler), false);
    for (std::size_t i = 0; i != Repeat; ++i) {
        if (Repeat > 1)
            std::cout << "Run: " << i + 1 << std::endl;
        zconst_file << schedule_name << ' ';
        zconst_file << alpha_config << ' ';
        smc_do(sampler, zconst_file, schedule_name, model_num, iter_num);
        // smc_do(sampler, zconst_file, schedule_name, CM, iter_num);
        zconst_file << std::endl;
    }
}

// modified from original to include model_num arguement
template <typename T, typename SD>
inline void smc_do (std::string &config, vsmc::Sampler<T> &sampler,
        std::ostream &zconst_file, std::size_t model_num)
{
    typedef move_smc<T, alpha_ess<T, cess01<T> >, SD> cess;
    typedef move_smc<T, alpha_prior    <T, 5>,    SD> prior5;
    typedef move_smc<T,alpha_fixed<T>, SD>user;

    if (config.compare("CESS")==0)
        for (std::size_t i = 0; i != CESSDrop.size(); ++i)
            smc_do<cess>(sampler, zconst_file, "CESS", CESSDrop[i], 0, model_num);

    if (config.compare("Prior5")==0){
        for (std::size_t i = 0; i != Prior5IterNum.size(); ++i)
            smc_do<prior5>(sampler, zconst_file, "Prior5",
                           Prior5IterNum[i], 0, model_num);
    }
    if (config.compare("User")==0){
        smc_do<user>(sampler, zconst_file, "User", &alpha_user_input, 0, model_num);
    }
    //#######################
    //Potentially interesting annealing schemes -- not tested
    // more schemes can be found in original vSMC.
    // typedef move_smc<T, alpha_mh<T>,              SD> mh;
    // typedef move_smc<T, alpha_ess<T, ess01<T> >,  SD> ess;
    // typedef move_smc<T, alpha_linear   <T>,       SD> linear;

    // if (config.count("mh_alpha") && config.count("mh_iter") &&
    //         MHAlpha.size() == MHIterNum.size())
    //     for (std::size_t i = 0; i != MHAlpha.size(); ++i)
    //         smc_do<mh>(sampler, zconst_file, "MH", MHAlpha[i], MHIterNum[i], model_num);
    //
    // if (config.count("ess"))
    //     for (std::size_t i = 0; i != ESSDrop.size(); ++i)
    //         smc_do<ess>(sampler, zconst_file, "ESS", ESSDrop[i], 0, model_num);
    //
    // if (config.count("linear"))
    //     for (std::size_t i = 0; i != LinearIterNum.size(); ++i)
    //         smc_do<linear>(sampler, zconst_file, "Linear",
    //                 LinearIterNum[i], 0, model_num);

}

#endif // VSMC_EXAMPLE_SMC_DO_HPP
