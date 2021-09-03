//============================================================================
// vSMC/example/pet/include/pet_init.hpp
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

#ifndef VSMC_EXAMPLE_PET_INIT_HPP
#define VSMC_EXAMPLE_PET_INIT_HPP
#include <RcppArmadillo.h>
using namespace Rcpp;

static const std::size_t InitCompNumBAR = 1;

class pet_init : public BASE_INIT<pet_state, pet_init>
{
    public :

    void pre_processor (vsmc::Particle<pet_state> &particle)
    {
        particle.value().alpha(0);
        particle.value().alpha_inc(0);
        particle.weight_set().set_equal_weight();
    }

    void initialize_param (vsmc::Particle<pet_state> &particle, void *info)
    {
        // std::cout<<"\n initialize_state param called !";
        if (particle.value().state(0, 0).comp_num() == 0)
            particle.value().comp_num(InitCompNum);
        if (info)
            particle.value().read_data(static_cast<const pet_info *>(info));
    }

    std::size_t initialize_state (vsmc::SingleParticle<pet_state> sp)
    {
        // std::cout<<"\n initialize_state called !";
        const std::size_t cn = sp.state(0).comp_num();
        // std::cout<<"\n comp_num DONE"<<cn;
        // double first = sp.particle().value().phi_lb0(0);
        // std::cout<<"\n fisrt DONE"<<first;
        for (std::size_t d = 0; d != cn; ++d) {
            double phi_lb0 = sp.particle().value().phi_lb0(d);
            double phi_ub0 = sp.particle().value().phi_ub0(d);
            // std::cout<<"\n phi_u and l done !";
            double theta_lb0 = d ?
                sp.state(0).theta(d - 1):
                sp.particle().value().theta_lb0(d);
            double theta_ub0 = sp.particle().value().theta_ub0(d);

            // vsmc::cxx11::uniform_real_distribution<> rphi(
            //         phi_lb0, phi_ub0);
            // vsmc::cxx11::uniform_real_distribution<> rtheta(
            //         theta_lb0, theta_ub0);
            // sp.state(0).phi(d) = rphi(sp.rng());
            // sp.state(0).theta(d) = rtheta(sp.rng());

            Rcpp::NumericVector proposal_phi = Rcpp::runif(1,phi_lb0,phi_ub0);
            Rcpp::NumericVector proposal_theta =Rcpp::runif(1,theta_lb0,theta_ub0);
            // std::cout<<"\n proposal_phi:"<<proposal_phi[0];
            // stop("STOPPED");
            sp.state(0).phi(d) = proposal_phi[0];
            sp.state(0).theta(d) = proposal_theta[0];

        }
        // std::cout<<"\n loop done !";
        double lambda_a0 = sp.particle().value().lambda_a0();
        double lambda_b0 = sp.particle().value().lambda_b0();
        // vsmc::cxx11::gamma_distribution<> rlambda(lambda_a0, lambda_b0);

        Rcpp::NumericVector proposal_lambda = Rcpp::rgamma(1,lambda_a0,lambda_b0);

        // sp.state(0).lambda() = rlambda(sp.rng());
        sp.state(0).lambda() = proposal_lambda[0];
        if (sp.state(0).lambda() < 1e-13)
            sp.state(0).lambda() = 1e-13;
        if (sp.particle().value().model() == StudentT) {
            double nu_a0 = sp.particle().value().nu_a0();
            double nu_b0 = sp.particle().value().nu_b0();
            // vsmc::cxx11::uniform_real_distribution<> rnu(nu_a0, nu_b0);
            Rcpp::NumericVector proposal_nu = Rcpp::runif(1,nu_a0,nu_b0);

            // sp.state(0).nu() = 1 / rnu(sp.rng());
            sp.state(0).nu() = 1 /proposal_nu[0];
        }
        sp.particle().value().log_target(sp.state(0));

        return 0;
    }
};

#endif // VSMC_EXAMPLE_PET_INIT_HPP
