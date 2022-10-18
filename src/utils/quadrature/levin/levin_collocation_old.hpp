#ifndef HEOM_TN_QUAD_LEVIN_BESSEL_COLLOCATION_HPP
#define HEOM_TN_QUAD_LEVIN_BESSEL_COLLOCATION_HPP

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

namespace heom_tn
{
namespace quad
{

template <typename T> 
class levin
{
    using RT = typename linalg::get_real_type<T>::type;
public:
    levin() : m_p(0) {}
    levin(size_t p, RT reg = 1e-10) : m_p(p), m_solver(2*p, true), m_zk(p), m_M(2*p, 2*p), m_f(2*p), m_c(2*p)
    {
        ASSERT(p > 1, "Levin collocation requires at least 2 points.");
        size_t N = m_p-1;
        for(size_t i=0; i<p; ++i)
        {
            m_zk(i) = std::cos(i*M_PI/N);
        }
    }

    levin(const levin& o) = default;
    levin(levin&& o) = default;
       
    levin& operator=(const levin& o) = default;
    levin& operator=(levin&& o) = default;

    void resize(size_t p)
    {
        ASSERT(p > 1, "Failed to resize Levin collocation object.  Levin collocation requires at least 2 points.");
        if(p != m_p)
        {
            m_p = p;
            m_zk.resize(p);
            m_M.resize(2*p, 2*p);
            m_f.resize(2*p);
            m_c.resize(2*p);
            m_solver.resize(2*p, true);

            size_t N = m_p-1;
            for(size_t i=0; i<p; ++i)
            {
                m_zk(i) = std::cos(i*M_PI/N);
            }
        }
    }

    //evaluate integral int_a^b f(x, args...)*J_nu(x*r) dx
    template <typename F, typename ... Args> 
    T operator()(F&& f, RT nu, RT r, RT a, RT b, Args&& ... args)
    {
        ASSERT(m_p > 1, "Cannot evaluate integral using levin collocation with 1 or fewer points.");
        CALL_AND_HANDLE(
            return eval_function
            (
                [](size_t i, size_t j, RT m, RT x, RT _r) -> RT 
                {
                    size_t ind = i*2+j;
                    switch(ind)
                    {
                        case(0):    
                            return m/x;
                        case(1):
                            return -_r;
                        case(2): 
                            return _r;
                        case(3):
                            return -(m+1)/x;
                    }
                },
                [](RT _nu, RT _r, RT x){return std::cyl_bessel_j(_nu, _r*x);},
                std::forward<F>(f),
                nu,
                r,
                a, 
                b, 
                std::forward<Args>(args)...
            ),
            "Failed to apply Levin collocation."
        ); 
    }       
    
    size_t size() const
    {
        return m_p;
    }
 
protected:
    template <typename F, typename A, typename W, typename ... Args> 
    T eval_function(A&& Af, W&& wf, F&& f, RT nu, RT r, RT a, RT b, Args&& ... args)
    {
        if(a > b)
        {
            RT c = a;
            a = b;
            b = c;
        }
        RT d = (b-a)/2.0;
        RT e = (b+a)/2.0;
        size_t N = m_p-1;

        //first evaluate the function at the collocation points
        for(size_t i=0; i<m_p; ++i)
        {
            RT x = m_zk(i)*d + e;
            m_f(i) = f(x, std::forward<Args>(args)...);         
            m_f(i+m_p) = 0.0;
        }

        //now construct the matrix that needs to be inverted.
        for(size_t i=0; i<m_p; ++i)
        {
            RT x = m_zk(i)*d + e;

            T A11 = Af(0, 0, nu, x, r);
            T A12 = Af(0, 1, nu, x, r);
            T A21 = Af(1, 0, nu, x, r);
            T A22 = Af(1, 1, nu, x, r);
            for(size_t k=0; k<m_p; ++k)
            {
                RT uij = cos(k*i*M_PI/N);
                RT upij = 0.0;
                if(k == 1)
                {
                    upij = 2*k/(b-a);
                }
                else if( k > 1)
                {
                    if(i == 0){             upij = 2*k*k/(b-a);}
                    else if (i+1 == m_p){   upij = 2*k*k/(b-a) * (k % 2 == 0 ? -1 : 1);}
                    else{                   upij = 2*k/(b-a)* sin(k*M_PI*i/N)/sin(M_PI*i/N);}
                }

                m_M(i, k) =  upij + A11*uij;
                m_M(i, k+m_p) =  A21*uij;
                m_M(i+m_p, k) =  A12*uij;
                m_M(i+m_p, k+m_p) = upij + A22*uij;
            }
        }

        //now solve the linear system of equations to find the coefficient matrix 
        CALL_AND_HANDLE(m_solver(m_M, m_f, m_c), "Failed to solve linear system.");

        std::array<T,2> Fija;   Fija[0] = 0.0;  Fija[1] = 0.0;
        std::array<T,2> Fijb;   Fijb[0] = 0.0;  Fijb[1] = 0.0;
        for(size_t i=0; i<2; ++i)
        {
            for(size_t k=0; k<m_p; ++k)
            {
                Fija[i] += (k % 2 == 0 ? 1.0 : -1.0) * m_c(k + i*m_p);
                Fijb[i] += m_c(k + i*m_p);
            }
        }

        return (Fijb[0] * wf(nu, r, b) + Fijb[1] * wf(nu+1, r, b)) - (Fija[0] * wf(nu, r, a) + Fija[1] * wf(nu+1, r, a));
    }


    size_t m_p;
    linalg::linear_solver<linalg::matrix<T> > m_solver;
    
    //store the chebyshev lobatto quadrature points 
    linalg::vector<RT> m_zk;

    linalg::matrix<T> m_M;
    linalg::vector<T> m_f;
    linalg::vector<T> m_c;
};

} //namespace quad
} //namespace heom_tn

#endif  //HEOM_TN_QUAD_LEVIN_BESSEL_COLLOCATION_HPP//

