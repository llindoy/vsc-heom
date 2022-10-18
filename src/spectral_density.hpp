#ifndef SUPER_OHMIC_SPECTRAL_DENSITY_HPP
#define SUPER_OHMIC_SPECTRAL_DENSITY_HPP

#include <linalg/linalg.hpp>
#include <vector>
#include <cmath>
#include <iomanip>

#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include "quadrature/adaptive_integrate.hpp"
#include "quadrature/gaussian_quadrature/gauss_legendre_quadrature.hpp"



template <typename T, typename Dens> 
inline T generate_frequencies(std::vector<T>& wk, Dens&& rho, T wmax, T tol = 1e-12, size_t max_iter = 1000)
{
    heom_tn::quad::gauss::legendre<T> leg(100);

    //determine the normalisation factor for the density so that \int_0^\omega_{max} \rho(\omega) \mathrm{d}\omega = N
    T c = heom_tn::quad::adaptive_integrate<T>([&](T w){return rho(w);}, leg, static_cast<T>(0.0), wmax, 1e-14, 0.0, true, 1e-14);
    size_t N = wk.size();
    T a1 = N/c;
    
    T w = wmax;
    for(size_t i=0; i < N; ++i)
    {
        size_t j = N-i;
        size_t ind = j-1;
    
        T delta_w = 10;
        T delta_f = 10;
    
        T df = rho(w)*a1;
        T f = heom_tn::quad::adaptive_integrate<T>([&](T _w){return rho(_w)*a1;}, leg, static_cast<T>(0.0), w, 1e-14, 0.0, true, 1e-14) - j;
    
        T _wmin = 0;
        T _wmax = w;
    
        delta_f = std::abs(f);
        T wtrial = w - f/df;
        if(wtrial < _wmin || wtrial > _wmax || std::isnan(wtrial))
        {
            T wm = w;
            if(f > 0){_wmax = w;    w = (w+_wmin)/2.0;}
            else{_wmin = w;         w = (w+_wmax)/2.0;}
            delta_w = std::abs(wm - w);
        }
        else
        {
            w = wtrial;
            delta_w = std::abs(f/df);
        }
        ASSERT_NUMERIC(!std::isnan(w), "Invalid frequency obtained in discretization of spectral_density.");
        size_t count = 0;
        
        while(delta_w > tol && delta_f > tol)  
        {
            T fx = heom_tn::quad::adaptive_integrate<T>([&](T _w){return rho(_w)*a1;}, leg, static_cast<T>(0.0), w, 1e-14, 0.0, true, 1e-14) - j;
            T dfx = rho(w)*a1; 
            delta_f = std::abs(fx);
            f = fx;
            df = dfx;
    
            wtrial = w - f/df;
            T wm = w;
    
            if(f > 0){_wmax = w;}
            else if(f < 0){_wmin = w;}
            if(wtrial < _wmin || wtrial > _wmax || std::isnan(wtrial))
            {
                if(f > 0){w = (w+_wmin)/2.0;}
                else{w = (w+_wmax)/2.0;}
            }
            else
            {
                if(f > 0){w = (w+_wmin)/2.0;}
                else{w = (w+_wmax)/2.0;}
                w = 0.1*w + 0.9*wtrial;
            }
            delta_w = std::abs(w - wm);
            ASSERT_NUMERIC(!std::isnan(w), "Invalid frequency obtained in discretization of spectral_density.");
            ++count;
            ASSERT(count < max_iter, "Newton's method failed to converge within the allowed number of iterators.");
        }
        wk[ind] = w;
    }
    return a1;
}

template <typename T> 
inline void debye_spectral_density(const T& Lambda, const T& wc, std::vector<T>& wk, std::vector<T>& gk)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    size_t N = wk.size();
    ASSERT(N == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");

    auto rho_lambda = [&Lambda, &wc](T w){return 2*Lambda*wc/(w*w+wc*wc);};
    auto J = [&Lambda, &wc](T w){return 2*Lambda*wc*w/(w*w+wc*wc);};
    T imax = N/(N+1.0)*M_PI/2.0;

    //for the debye spectral density we want to choose our maximum frequency so that we capture (N-1)/N*\lambda of the renormalisation energy
    T wmax = wc*std::tan(imax); 
    T a1 = generate_frequencies(wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        gk[i] =  sqrt(1.0/M_PI*J(wk[i])/(a1*rho_lambda(wk[i])));
    }
}

template <typename T> 
inline void debye_spectral_density(const T& Lambda, const T& wc, const T& beta, std::vector<T>& wk, std::vector<T>& gk,  T renorm = 1e-2)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    ASSERT(wk.size() == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    size_t N = wk.size()/2;


    auto J = [&Lambda, &wc, &beta](T w)
    {
        T x = std::abs(w);
        return (w > 0 ? 1.0 : -1.0)*Lambda*wc*x/(x*x+wc*wc) * (1.0 + 1.0/std::tanh(beta*w/2.0));
    };
    auto rho_lambda = [&Lambda, &wc, &beta, &renorm](T w)
    {
        return Lambda*wc/(w*w+wc*wc) * (1.0 + std::tanh(beta*w/2.0)/(renorm*renorm+std::tanh(beta*w/2.0)*std::tanh(beta*w/2.0)));
    };

    std::vector<T> _wk(N);

    T imax = N/(N+1.0)*M_PI/2.0;
    //for the debye spectral density we want to choose our maximum frequency so that we capture (N-1)/N*\lambda of the renormalisation energy.
    //The expression for the renormalization energy as a function of the upper bound is \propto atan(x/wc).  So we can analytically evaluate wmax
    T wmax = wc*std::tan(imax);
    T a1 = generate_frequencies(_wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        //wk[N-(i+1)] = -_wk[i];
        //wk[N+i] = _wk[i];

        //gk[N-(i+1)] = sqrt(2.0/M_PI*_wk[i]*Jm_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));
        //gk[N+i] =  sqrt(2.0/M_PI*_wk[i]*Jp_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));
        wk[2*i] = -_wk[i];
        wk[2*i+1] = _wk[i];

        gk[2*i] = sqrt(1.0/M_PI*J(wk[2*i])/(a1*rho_lambda(_wk[i])));
        gk[2*i+1] =  sqrt(1.0/M_PI*J(wk[2*i+1])/(a1*rho_lambda(_wk[i])));
    }
}


template <typename T> 
inline void brownian_oscillator_spectral_density(const T& Lambda, const T& Omega, const T& gamma, std::vector<T>& wk, std::vector<T>& gk)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    size_t N = wk.size();
    ASSERT(N == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");


    auto rho_lambda = [&Lambda, &Omega, &gamma](T w){return Lambda/2*gamma*Omega*Omega/((w*w-Omega*Omega)*(w*w-Omega*Omega)+gamma*gamma*w*w);};
    T wmax = Omega*2*std::log(N+1.0);//find_maximum_frequency(Lambda, std::forward<decltype(rho_lambda)>(rho_lambda), itol);
    T a1 = generate_frequencies(wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        gk[i] = sqrt(1.0/M_PI*wk[i]/a1);
    }
}

template <typename T> 
inline void brownian_oscillator_spectral_density(const T& Lambda, const T& Omega, const T& gamma, const T& beta, std::vector<T>& wk, std::vector<T>& gk, T renorm = 1e-2)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    ASSERT(wk.size() == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    size_t N = wk.size()/2;


    auto J = [&Lambda, &Omega, &gamma, &beta](T w)
    {
        T x = std::abs(w);
        return (w > 0 ? 1.0 : -1.0)*(Lambda/4*gamma*Omega*Omega*x/((x*x-Omega*Omega)*(x*x-Omega*Omega)+gamma*gamma*x*x)) * (1.0 + 1.0/std::tanh(beta*w/2.0));
    };
    auto rho_lambda = [&Lambda, &Omega, &gamma, &beta, &renorm](T w)
    {
        return (Lambda/2*gamma*Omega*Omega/((w*w-Omega*Omega)*(w*w-Omega*Omega)+gamma*gamma*w*w)) * 0.5 * (1.0 + std::tanh(beta*w/2.0)/(renorm*renorm+std::tanh(beta*w/2.0)*std::tanh(beta*w/2.0)));
    };

    std::vector<T> _wk(N);
    T wmax = Omega*2*std::log(N+1.0);//find_maximum_frequency(Lambda, std::forward<decltype(rho_lambda)>(rho_lambda), itol);
    T a1 = generate_frequencies(_wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        //wk[N-(i+1)] = -_wk[i];
        //wk[N+i] = _wk[i];

        //gk[N+i] =  sqrt(2.0/M_PI*_wk[i]*Jp_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));
        //gk[N-(i+1)] = sqrt(2.0/M_PI*_wk[i]*Jm_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));

        wk[2*i] = -_wk[i];
        wk[2*i+1] = _wk[i];

        gk[2*i] = sqrt(1.0/M_PI*J(wk[2*i])/(a1*rho_lambda(_wk[i])));
        gk[2*i+1] =  sqrt(1.0/M_PI*J(wk[2*i+1])/(a1*rho_lambda(_wk[i])));
    }
}


template <typename T> 
inline void ohmic_spectral_density(const T& alpha, const T& wc, std::vector<T>& wk, std::vector<T>& gk)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    size_t N = wk.size();
    ASSERT(N == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");


    auto rho_lambda = [&alpha, &wc](T w){return M_PI/2.0*alpha*exp(-w/wc);};
    T wmax = wc*std::log(N+1.0);//find_maximum_frequency(Lambda, std::forward<decltype(rho_lambda)>(rho_lambda), itol);
    T a1 = generate_frequencies(wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        gk[i] = sqrt(1.0/M_PI*wk[i]/a1);
    }
}

template <typename T> 
inline void ohmic_spectral_density(const T& alpha, const T& wc, const T& beta, std::vector<T>& wk, std::vector<T>& gk, T renorm = 1)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    ASSERT(wk.size() == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    size_t N = wk.size()/2;

    auto J = [&alpha, &wc, &beta](T w)
    {
        T x = std::abs(w);
        return M_PI/4.0*(w > 0 ? 1.0 : -1.0)*alpha*x*exp(-x/wc)* (1.0 + 1.0/std::tanh(beta*w/2.0));
    };
    auto rho_lambda = [&alpha, &wc, &beta, &renorm](T w)
    {
        return M_PI/2.0*alpha*exp(-w/wc) * 0.5 * (1.0 + std::tanh(beta*w/2.0)/(renorm*renorm+std::tanh(beta*w/2.0)*std::tanh(beta*w/2.0)));
    };

    std::vector<T> _wk(N);
    T wmax = wc*std::log(N+1.0);//find_maximum_frequency(Lambda, std::forward<decltype(rho_lambda)>(rho_lambda), itol);
    T a1 = generate_frequencies(_wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        //wk[N-(i+1)] = -_wk[i];
        //wk[N+i] = _wk[i];

        //gk[N+i] =  sqrt(2.0/M_PI*_wk[i]*Jp_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));
        //gk[N-(i+1)] = sqrt(2.0/M_PI*_wk[i]*Jm_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));

        wk[2*i] = -_wk[i];
        wk[2*i+1] = _wk[i];

        gk[2*i] = sqrt(1.0/M_PI*J(wk[2*i])/(a1*rho_lambda(_wk[i])));
        gk[2*i+1] =  sqrt(1.0/M_PI*J(wk[2*i+1])/(a1*rho_lambda(_wk[i])));
    }
}


template <typename T, size_t D> struct spec_dense;

template <typename T> struct spec_dense<T, 1>
{
    static T eval(T x)
    {
        return cos(x);
    }

    static T distance(const std::array<T, 1>& a, const std::array<T, 1>& b)
    {
        return std::abs(a[0]-b[0]);
    }

    static T gamma(T wc, T wm, size_t N)
    {
        return static_cast<T>(N)/(wc - wc*exp(-wm/wc));
    }

    static T fopt(size_t j, T g, T wc, T w)
    {
        return wc - wc*exp(-w/wc) - static_cast<T>(j)/g;
    }

    static T fderiv(T wc, T w)
    {
        return exp(-w/wc);
    }
    static T fscale(T)
    {
        return 1.0;
    }
};

/*  
template <typename T> struct spec_dense<T, 2>
{
    static T eval(T x)
    {
        return std::cyl_bessel_j(static_cast<T>(0.0), x);
    }

    static T distance(const std::array<T, 2>& a, const std::array<T, 2>& b)
    {
        std::array<T, 2> c;
        for(size_t i=0; i < 2; ++i){c[i] = a[i]-b[i];}
        return sqrt(c[0]*c[0]+c[1]*c[1]);
    }

    static T gamma(T wc, T wm, size_t N)
    {
        
        return static_cast<T>(N)/(wc - (wc + wm)*exp(-wm/wc));
    }

    static T fopt(size_t j, T g, T wc, T w)
    {
        return wc - (wc+w)*exp(-w/wc) - static_cast<T>(j)/g;
    }

    static T fderiv(T wc, T w)
    {
        return w*exp(-w/wc)/wc;
    }
};
*/
template <typename T> struct spec_dense<T, 3> 
{
    static T eval(T x)
    {
        if(x < 1e-5)
        {
            T x2 = x*x;
            T x4 = x2*x2;
            T x6 = x4*x2;
            T x8 = x4*x4;
            return 1 - x2/6.0 + x4/120 - x6/5040 + x8/362880;
        }
        else
        {
            return sin(x)/x;
        }
    }

    static T distance(const std::array<T, 3>& a, const std::array<T, 3>& b)
    {
        std::array<T, 3> c;
        for(size_t i=0; i < 3; ++i){c[i] = a[i]-b[i];}
        return sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
    }

    static T gamma(T wc, T wm, size_t N)
    {
        return (N/(2*wc*wc - (2*wc*wc + 2*wc*wm + wm*wm)*std::exp(-wm/wc)));
    }

    static T fopt(size_t j, T g, T wc, T w)
    {
        return (2.0*wc*wc - static_cast<T>(j)/(g)) - (2*wc*wc+2*wc*w+w*w)*exp(-w/wc);
    }

    static T fderiv(T wc, T w)
    {
        return w*w/wc*exp(-w/wc);
    }

    static T fscale(T wc)
    {
        return wc;
    }
};

template <typename T> 
T Jzzzz(T x)
{
    T x2 = x*x;
    T x4 = x2*x2;
    if(x < 1e-5){return 1 - 5*x2/14.0 + 5*x2*x2/216 - x4*x2/1584 + x4*x4/104832;}
    else{return 5*(4*x*(x2-6)*cos(x) + (x4 - 12*x2+24)*sin(x))/(x4*x);}
}


template <typename T>
inline void cholesky_decompose(linalg::matrix<T>& A, linalg::matrix<T>& L)
{
    ASSERT(A.shape(0) == A.shape(1) && A.shape(0) == L.shape(0) && A.shape(1) == L.shape(1), "Invalid matrices for cholesky decompose.");
    L.fill_zeros();

    auto hm = hermitian_view(A);
    linalg::eigensolver<decltype(hm)> solver(hm);

    linalg::diagonal_matrix<T> vals;
    linalg::matrix<T> vecs;
    solver(hm, vals, vecs);
    for(size_t i=0; i < vals.size(); ++i)
    {
        if(std::abs(vals(i, i)/vals(vals.size()-1, vals.size()-1)) < 1e-12){vals(i, i) = 0.0;}
        else{vals(i, i) = sqrt(vals(i, i));}
    }

    A = vecs*vals;
    L = A*adjoint(vecs);
}


template <typename T, size_t D> struct correlation;
template <typename T> struct correlation<T, 1>
{
public:
    static T f(const T& w, const std::array<T, 1>& a, const std::array<T, 1>& b)
    {
        return std::cos(w*std::abs(a[0]-b[0]));
    }
};

/*  
template <typename T, size_t D> 
inline void correlated_spectral_density(const std::vector<T>& alpha, const std::vector<std::array<T, D>>& x, const T& wc, std::vector<T>& wk, std::vector<linalg::matrix<T>>& gk)
{
    size_t nb = alpha.size();
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");

    size_t N = wk.size();
    ASSERT(N == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    ASSERT(nb > 0 , "Invalid number of spins.");
    ASSERT(nb == x.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    
    for(size_t i = 0; i < N; ++i)
    {
        ASSERT(gk[i].shape(0) == nb, "Failed to compute the discretised bath spectral density.  Invalid sizes.");
        ASSERT(gk[i].shape(1) == nb, "Failed to compute the discretised bath spectral density.  Invalid sizes.");
    }

    linalg::matrix<T> Jr(nb, nb);

    auto rho_lambda = [&alpha, &wc](T w){return M_PI/2.0*exp(-w/wc);};
    T wmax = wc*std::log(N+1.0);//find_maximum_frequency(Lambda, std::forward<decltype(rho_lambda)>(rho_lambda), itol);
    T a1 = generate_frequencies(wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        for(size_t bi=0; bi < nb; ++bi)
        {
            Jr(bi, bi) = 1.0/M_PI*wk[i]*alpha[bi]/(a1);
            for(size_t bj=0; bj < bi; ++bj)
            {
                Jr(bi, bj) = sqrt(Jr(bi, bi)*Jr(bj, bj)) * correlation<T, D>::eval(wk[i], x[bi], x[bj]);
                Jr(bj, bi) = Jr(bi, bj);
            }
        }
        cholesky_decompose(Jr, gk[i]);
    }
}

template <typename T, size_t D> 
inline void correlated_spectral_density(const std::vector<T>& alpha, const std::vector<std::array<T, D>>& x, const T& wc, const T& beta, std::vector<T>& wk, std::vector<linalg::matrix<T>>& gk, T renorm = 1)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    ASSERT(wk.size() == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");

    size_t N = wk.size()/2;
    ASSERT(wk.size() == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    ASSERT(nb > 0 , "Invalid number of spins.");
    ASSERT(nb == x.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    
    for(size_t i = 0; i < N; ++i)
    {
        ASSERT(gk[i].shape(0) == nb, "Failed to compute the discretised bath spectral density.  Invalid sizes.");
        ASSERT(gk[i].shape(1) == nb, "Failed to compute the discretised bath spectral density.  Invalid sizes.");
    }

    auto J = [&alpha, &wc, &beta](T w)
    {
        T x = std::abs(w);
        return M_PI/4.0*(w > 0 ? 1.0 : -1.0)*x*exp(-x/wc)* (1.0 + 1.0/std::tanh(beta*w/2.0));
    };
    auto rho_lambda = [&alpha, &wc, &beta, &renorm](T w)
    {
        return M_PI/2.0*exp(-w/wc) * 0.5 * (1.0 + std::tanh(beta*w/2.0)/(renorm*renorm+std::tanh(beta*w/2.0)*std::tanh(beta*w/2.0)));
    };

    std::vector<T> _wk(N);
    T wmax = wc*std::log(N+1.0);//find_maximum_frequency(Lambda, std::forward<decltype(rho_lambda)>(rho_lambda), itol);
    T a1 = generate_frequencies(_wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {

        for(size_t bi=0; bi < nb; ++bi)
        {
            Jr(bi, bi) = 1.0/M_PI*alpha[bi]*J(wk[2*i])/(a1*rho_lambda(_wk[i]));
            for(size_t bj=0; bj < bi; ++bj)
            {
                Jr(bi, bj) = sqrt(Jr(bi, bi)*Jr(bj, bj)) * correlation<T, D>::eval(wk[i], x[bi], x[bj]);
                Jr(bj, bi) = Jr(bi, bj);
            }
        }
        cholesky_decompose(Jr, gk[2*i]);


        for(size_t bi=0; bi < nb; ++bi)
        {
            Jr(bi, bi) = 1.0/M_PI*alpha[bi]*J(wk[2*i+1])/(a1*rho_lambda(_wk[i]));
            for(size_t bj=0; bj < bi; ++bj)
            {
                Jr(bi, bj) = sqrt(Jr(bi, bi)*Jr(bj, bj)) * correlation<T, D>::eval(wk[i], x[bi], x[bj]);
                Jr(bj, bi) = Jr(bi, bj);
            }
        }
        cholesky_decompose(Jr, gk[2*i+1]);
    }
}*/

template <typename T> 
T super_ohmic_max_freq(const T& wc, size_t N, T tol = 1e-12)
{
    T _xmin = 0;
    
    T x = 10;
    T _xmax = 20;
    while(std::exp(-_xmax)*(_xmax+1) > 1.0/(N+1.0))
    {
        _xmax += 20;
    }

    T delta_x = 1e6;
    T delta_f = 1e6;
    T f, df, xtrial;
    size_t count = 0;
    while(delta_x > tol && delta_f > tol)  
    {
        T fx = std::exp(-x)*(x+1) - 1.0/(N+1.0);
        T dfx = -x*std::exp(-x); 
        delta_f = std::abs(fx);
        f = fx;
        df = dfx;
    
        xtrial = x - f/df;
        T xm = x;
    
        if(f > 0){_xmax = x;}
        else if(f < 0){_xmin = x;}
        if(xtrial < _xmin || xtrial > _xmax || std::isnan(xtrial))
        {
            if(f > 0){x = (x+_xmin)/2.0;}
            else{x = (x+_xmax)/2.0;}
        }
        else
        {
            if(f > 0){x = (x+_xmin)/2.0;}
            else{x = (x+_xmax)/2.0;}
            x = 0.1*x + 0.9*xtrial;
        }
        delta_x = std::abs(x - xm);
        ASSERT_NUMERIC(!std::isnan(x), "Invalid frequency obtained in discretization of spectral_density.");
        ++count;
        ASSERT(count < 1000, "Newton's method failed to converge within the allowed number of iterators.");
    }
    return x*wc;
}


template <typename T> 
inline void super_ohmic_spectral_density(const T& alpha, const T& wc, std::vector<T>& wk, std::vector<T>& gk)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    size_t N = wk.size();
    ASSERT(N == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");


    auto rho_lambda = [&alpha, &wc](T w){return M_PI/2.0*alpha/wc*w*exp(-w/wc);};
    T wmax = super_ohmic_max_freq(wc, N);
    T a1 = generate_frequencies(wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        gk[i] = sqrt(1.0/M_PI*wk[i]/a1);
    }
}

template <typename T> 
inline void super_ohmic_spectral_density(const T& alpha, const T& wc, const T& beta, std::vector<T>& wk, std::vector<T>& gk, T renorm = 1e-2)
{
    ASSERT(wk.size() > 0, "Failed to compute the discretised bath spectral density.  Invalid size");
    ASSERT(wk.size() == gk.size(), "Failed to compute the discretised bath spectral density.  The two input arrays are not the same size.");
    size_t N = wk.size()/2;


    auto J = [&alpha, &wc, &beta](T w)
    {
        T x = std::abs(w);
        return M_PI/4.0*(w > 0 ? 1.0 : -1.0)*alpha*x*x*exp(-x/wc)* (1.0 + 1.0/std::tanh(beta*w/2.0));
    };
    auto rho_lambda = [&alpha, &wc, &beta, &renorm](T w)
    {
        if(std::abs(beta*w/2) > 1e-4)
        {
            return M_PI/4.0*alpha*w/wc*exp(-w/wc) * (1.0 + 1.0/std::tanh(beta*w/2.0));
        }
        else
        {
            T wbose = w + 2.0/beta;
            for(size_t i = 1; i < 1000; ++i)
            {
                wbose += beta*w*w/(M_PI*M_PI*i*i+beta*beta*w*w/4);
            }
            return M_PI/4.0*alpha/wc*exp(-w/wc) * wbose;
        }
    };

    std::vector<T> _wk(N);
    T wmax = super_ohmic_max_freq(wc, N);
    T a1 = generate_frequencies(_wk, std::forward<decltype(rho_lambda)>(rho_lambda), wmax);

    for(size_t i = 0; i < N; ++i)
    {
        //wk[N-(i+1)] = -_wk[i];
        //wk[N+i] = _wk[i];

        //gk[N-(i+1)] = sqrt(2.0/M_PI*_wk[i]*Jm_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));
        //gk[N+i] =  sqrt(2.0/M_PI*_wk[i]*Jp_lambda(_wk[i])/(a1*rho_lambda(_wk[i])));
        wk[2*i] = -_wk[i];
        wk[2*i+1] = _wk[i];

        gk[2*i] = sqrt(1.0/M_PI*J(wk[2*i])/(a1*rho_lambda(_wk[i])));
        gk[2*i+1] =  sqrt(1.0/M_PI*J(wk[2*i+1])/(a1*rho_lambda(_wk[i])));
    }
}

#endif

