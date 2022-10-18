#ifndef SUPER_OHMIC_SPECTRAL_DENSITY_HPP
#define SUPER_OHMIC_SPECTRAL_DENSITY_HPP

#include <linalg/linalg.hpp>
#include <vector>
#include <cmath>
#include <iomanip>

#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include "quadrature/adaptive_integrate.hpp"
#include "quadrature/gaussian_quadrature/gauss_legendre_quadrature.hpp"



template <typename T> 
inline void debye_correlation_function(const T& Lambda, const T& wc, const T& beta, T& delta, std::vector<linalg::complex<T>>& ck, std::vector<linalg::complex<T>>& nuk)
{
    ASSERT(ck.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    size_t N = ck.size();
    ASSERT(N == nuk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");

    auto Ck = [&Lambda, &wc, &beta](T w){return 4*Lambda*wc/beta/(w*w-wc*wc);};

    nuk[0] = linalg::complex<T>(wc, 0);
    ck[0] = Lambda*wc*linalg::complex<T>(1.0/std::tan(beta*wc/2.0),-1);
    for(size_t k=1; k < N; ++k)
    {
        T nu = 2*M_PI*k/beta; 
        nuk[k] = linalg::complex<T>(nu, 0);
        ck[k] = Ck(nu)*nuk[k]; 
    }
    delta = 2.0*Lambda/(beta*wc);
    for(size_t k=0; k < N-1; ++k)
    {
        delta -= std::real(ck[k])/std::real(nuk[k]);
    }
}



template <typename T> 
inline void brownian_oscillator_correlation_function(const T& Lambda, const T& Omega, const T& gamma, const T& beta, T& delta, std::vector<linalg::complex<T>>& ck, std::vector<linalg::complex<T>>& nuk)
{
    ASSERT(ck.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    size_t N = ck.size();
    ASSERT(N == nuk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");

    ASSERT(std::abs(2*Omega/gamma - 1) > 1e-12, "We do not treat the critically damped case here.");

    auto Ck = [&Lambda, &Omega, &gamma](T w){return 2*Lambda*gamma*Omega*Omega/(gamma*gamma*w*w-(w*w+gamma*gamma)*(w*w+gamma*gamma));};

    //underdamped
    if(2*Omega > gamma)
    {
        nuk[0] = linalg::complex<T>(gamma/2, sqrt(Omega*Omega - gamma*gamma/4));
        nuk[1] = linalg::complex<T>(gamma/2,-sqrt(Omega*Omega - gamma*gamma/4));

        //ck[0] = 
        //ck[1] =
    }
    //overdamped
    else
    {
        nuk[0] = gamma/2 + sqrt(gamma*gamma/4-Omega*Omega);
        nuk[1] = gamma/2 - sqrt(gamma*gamma/4-Omega*Omega);

        //ck[0] = 
        //ck[1] =
    }
    for(size_t k=2; k < N; ++k)
    {
        size_t ik = k-1;
        T nu = 2*M_PI*ik/beta; 
        nuk[k] = nu;
        ck[k] = 2/beta*Ck(nu)*nuk[k]; 
    }
    for(size_t k=N; k < 1e6; ++k)
    {
        size_t ik = k-1;
        T nu = 2*M_PI*ik/beta; 
        delta += Ck(nu)*2/beta;
    }

}

#endif

