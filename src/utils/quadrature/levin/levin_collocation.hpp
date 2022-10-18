#ifndef QUAD_LEVIN_COLLOCATION_HPP
#define QUAD_LEVIN_COLLOCATION_HPP

//A Levin collocation routine for evaluating oscillatory integrals involving functions multiplied by bessel functions.
//This approach takes the problem of numerically integrating the oscillatory function and transforms it onto an inverse problem.  
//Additionally this approach is expected to become more efficient as the frequency of oscillation increases.

#include <cstdint>
#include <cmath>
#include <array>
#include <limits>
#include <linalg/linalg.hpp>
#include <linalg/special_functions/linear_solver.hpp>
#include <linalg/decompositions/singular_value_decomposition/singular_value_decomposition.hpp>

namespace quad
{

template <typename Impl> class levin_traits;

template <class Impl>
class levin_base
{
public:
    using result_type = typename levin_traits<Impl>::result_type;
    using value_type = typename levin_traits<Impl>::value_type;
    using real_type = typename levin_traits<Impl>::real_type;
public:
    levin_base() : m_p(0) {}
    levin_base(size_t p) : m_p(0)
    {       
        CALL_AND_HANDLE(resize(p), "Failed to construct levin base object.");
    }

    levin_base(const levin_base& o) = default;
    levin_base(levin_base&& o) = default;
       
    levin_base& operator=(const levin_base& o) = default;
    levin_base& operator=(levin_base&& o) = default;
    
    void resize(size_t p)
    {
        try
        {
            if(m_p != p)
            {
                ASSERT(p > 1, "Must have at least 2 collocation points to employ levin collocation.");
                m_p = p;
                m_n = Impl::nfuncs()*p;
                CALL_AND_HANDLE(m_solver.resize(m_n, true), "Failed to resize linear solver.");
                CALL_AND_HANDLE(m_zk.resize(m_p), "Failed to resize collocation point array.");
                CALL_AND_HANDLE(m_M.resize(m_n, m_n), "Failed to resize M matrix.");
                CALL_AND_HANDLE(m_f.resize(m_n), "Failed to resize array storing value of function at collocation points.");
                CALL_AND_HANDLE(m_c.resize(m_n), "Failed to resize coefficient array.");

                //evaluate the chebyshev collocation points
                size_t N = m_p-1;   for(size_t i=0; i<p; ++i){m_zk(i) = std::cos(i*M_PI/N);}
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to resize levin_base object.");
        }
    }

    result_type oscillatory_function(real_type x){return static_cast<Impl*>(this)->function(0, x);}

    template <typename F, typename ... Args> 
    result_type operator()(F&& f, real_type r, real_type a, real_type b, Args&& ... args)
    {
        try
        {
            ASSERT(m_p > 1, "Cannot evaluate integral using levin collocation with 1 or fewer points.");

            if(a > b)
            {
                real_type c = a;
                a = b;
                b = c;
            }
            real_type d = (b-a)/2.0;
            real_type e = (b+a)/2.0;
            size_t N = m_p-1;

            //first evaluate the function at the collocation points
            m_f.fill_zeros();
            m_M.fill_zeros();
            for(size_t i=0; i<m_p; ++i)
            {
                real_type x = m_zk(i)*d + e;
                m_f(i) = f(x, std::forward<Args>(args)...);         
            }

            std::array<std::array<result_type, Impl::nfuncs()>, Impl::nfuncs()> A;

            //now construct the matrix that needs to be inverted.
            for(size_t i=0; i<m_p; ++i)
            {
                real_type x = m_zk(i)*d + e;

                for(size_t j=0; j<Impl::nfuncs(); ++j)
                {
                    for(size_t k=0; k<Impl::nfuncs(); ++k)
                    {
                        A[j][k] = static_cast<Impl*>(this)->differential_matrix(j, k, x, r);
                    }
                }

                for(size_t k=0; k<m_p; ++k)
                {
                    real_type uij = cos(k*i*M_PI/N);
                    real_type upij = 0.0;
                    if(k == 1){upij = 2*k/(b-a);}
                    else if( k > 1)
                    {
                        if(i == 0){             upij = 2*k*k/(b-a);}
                        else if (i+1 == m_p){   upij = 2*k*k/(b-a) * (k % 2 == 0 ? -1 : 1);}
                        else{                   upij = 2*k/(b-a)* sin(k*M_PI*i/N)/sin(M_PI*i/N);}
                    }

                    for(size_t m=0; m<Impl::nfuncs(); ++m)
                    {
                        m_M(i+m*m_p, k+m*m_p) = upij;
                        for(size_t n=0; n<Impl::nfuncs(); ++n)
                        {
                            m_M(i+m*m_p, k+n*m_p) += A[n][m]*uij;
                        }
                    }
                }
            }

            //now solve the linear system of equations to find the coefficient matrix 
            CALL_AND_HANDLE(m_solver(m_M, m_f, m_c), "Failed to solve linear system.");

            std::array<result_type,Impl::nfuncs()> Fija;  
            std::array<result_type,Impl::nfuncs()> Fijb;  

            for(size_t i=0; i<Impl::nfuncs(); ++i)
            {
                Fija[i] = 0.0;  Fijb[i] = 0.0;
                for(size_t k=0; k<m_p; ++k)
                {
                    Fija[i] += (k % 2 == 0 ? 1.0 : -1.0) * m_c(k + i*m_p);
                    Fijb[i] += m_c(k + i*m_p);
                }
            }

            result_type ret(0.0);
            
            for(size_t i=0; i<Impl::nfuncs(); ++i)
            {
                ret += Fijb[i] * static_cast<Impl*>(this)->function(i, r*b) - Fija[i] * static_cast<Impl*>(this)->function(i, r*a);
            }
            return ret;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to evaluate oscillatory integral using levin collocation.");
        }
    }

protected:

    size_t m_p, m_n;
    linalg::linear_solver<linalg::matrix<result_type> > m_solver;
    
    //store the chebyshev lobatto quadrature points 
    linalg::vector<real_type> m_zk;

    linalg::matrix<result_type> m_M;
    linalg::vector<result_type> m_f;
    linalg::vector<result_type> m_c;
};


template <typename T>
class levin_bessel : public levin_base<levin_bessel<T> >
{
public:
    using result_type = T;
    using real_type = typename linalg::get_real_type<T>::type;
    using value_type = T;

    using base_type = levin_base<levin_bessel<T>>;
    using levin_base<levin_bessel<T>>::operator=;
public:
    levin_bessel() : base_type(), m_nu(0) {}
    levin_bessel(size_t p) 
    try : base_type(p), m_nu(0) {}
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to construct levin_bessel object.");
    }

    levin_bessel(const levin_bessel& o) = default;
    levin_bessel(levin_bessel&& o) = default;

    static constexpr size_t nfuncs(){return 2;}
    value_type differential_matrix(size_t i, size_t j, real_type x, real_type r)
    {
        ASSERT(i < nfuncs() && j < nfuncs(), "Invalid input to differential matrix function.");
        size_t ind = i*nfuncs()+j;
        switch(ind)
        {
            case(0):    
                return m_nu/x;
            case(1):
                return -r;
            case(2): 
                return r;
            case(3):
                return -(m_nu+1)/x;
        };
        return 0;
    }
    
    result_type function(size_t i, real_type x)
    {
        ASSERT(i < nfuncs(), "Invalid input to function.");
        return std::cyl_bessel_j(m_nu+i, x);
    }

    real_type& nu(){return m_nu;}
    const real_type& nu() const{return m_nu;}
protected:
    real_type m_nu;
};

template <typename T>
struct levin_traits<levin_bessel<T>>
{
    using result_type = T;
    using real_type = typename linalg::get_real_type<T>::type;
    using value_type = T;
};

template <typename T>
class levin_fourier : public levin_base<levin_fourier<T> >
{
public:
    using real_type = typename linalg::get_real_type<T>::type;
    using result_type = std::complex<real_type>;
    using value_type = T;

    using base_type = levin_base<levin_fourier<T>>;
    using base_type::operator=;
    using base_type::levin_base;
public:

    static constexpr size_t nfuncs(){return 1;}
    result_type differential_matrix(size_t i, size_t j, real_type /*x*/, real_type r)
    {
        ASSERT(i < nfuncs() && j < nfuncs(), "Invalid input to differential matrix function.");
        return result_type(0, r);
    }
    
    result_type function(size_t i, real_type x)
    {
        ASSERT(i < nfuncs(), "Invalid input to function.");
        return exp(result_type(0, x));
    }
};

template <typename T>
struct levin_traits<levin_fourier<T>>
{
    using real_type = typename linalg::get_real_type<T>::type;
    using result_type = std::complex<real_type>;
    using value_type = T;
};

template <typename T>
class levin_sine : public levin_base<levin_sine<T> >
{
public:
    using real_type = typename linalg::get_real_type<T>::type;
    using result_type = T;
    using value_type = T;

    using base_type = levin_base<levin_sine<T>>;
    using base_type::operator=;
    using base_type::levin_base;
public:

    static constexpr size_t nfuncs(){return 2;}
    result_type differential_matrix(size_t i, size_t j, real_type /*x*/, real_type r)
    {
        ASSERT(i < nfuncs() && j < nfuncs(), "Invalid input to differential matrix function.");
        if(i == j){return 0;}
        else{return (i < j ? 1.0 : -1.0) * r;}
    }
    
    result_type function(size_t i, real_type x)
    {
        ASSERT(i < nfuncs(), "Invalid input to function.");
        if(i == 0){return sin(x);}
        else{return cos(x);}
    }
};

template <typename T>
struct levin_traits<levin_sine<T>>
{
    using result_type = T;
    using real_type = typename linalg::get_real_type<T>::type;
    using value_type = T;
};

template <typename T>
class levin_cosine : public levin_base<levin_cosine<T> >
{
public:
    using real_type = typename linalg::get_real_type<T>::type;
    using result_type = T;
    using value_type = T;

    using base_type = levin_base<levin_cosine<T>>;
    using base_type::operator=;
    using base_type::levin_base;
public:

    static constexpr size_t nfuncs(){return 2;}
    result_type differential_matrix(size_t i, size_t j, real_type /*x*/, real_type r)
    {
        ASSERT(i < nfuncs() && j < nfuncs(), "Invalid input to differential matrix function.");
        if(i == j){return 0;}
        else{return (i > j ? 1.0 : -1.0) * r;}
    }
    
    result_type function(size_t i, real_type x)
    {
        ASSERT(i < nfuncs(), "Invalid input to function.");
        if(i == 0){return cos(x);}
        else{return sin(x);}
    }
};

template <typename T>
struct levin_traits<levin_cosine<T>>
{
    using result_type = T;
    using real_type = typename linalg::get_real_type<T>::type;
    using value_type = T;
};

} //namespace quad

#endif  //HEOM_TN_QUAD_LEVIN_BESSEL_COLLOCATION_HPP//

