#ifndef GAUSS_LEGENDRE_QUADRATURE_HPP
#define GAUSS_LEGENDRE_QUADRATURE_HPP

#include <cstdint>
#include <linalg/linalg.hpp>
#include <linalg/decompositions/eigensolvers/eigensolver.hpp>

namespace quad{
namespace gauss{
template<typename T, typename backend>
class legendre_quad;

template <typename T>
class legendre_quad<T, linalg::blas_backend>
{
    using RT = typename linalg::get_real_type<T>::type;
public:
    template <typename F> 
    static auto impl(F&& f, const RT& a, const RT& b, const linalg::vector<RT, linalg::blas_backend>& w, const linalg::vector<RT, linalg::blas_backend>& x) -> decltype(f(b))
    {
        decltype(f(b)) ret = 0.0;
        for(size_t i=0; i<w.size(); ++i)
        {
            RT xp = (b-a)*(x(i)+1.0)/2.0 + a;
            ret += f(xp)*w(i);
        }
        ret *= (b-a)/2.0;
        return ret;
    }
};

//constructs and applies a gaussian quadrature rule of a given order for a specific orthogonal polynomial
template <typename T, typename backend = linalg::blas_backend>
class legendre
{
    using RT = typename linalg::get_real_type<T>::type;
public:
    legendre() {} 
    legendre(size_t N) : m_w(N), m_x(N)
    {
        construct_quadrature_rule();
    }
    ~legendre(){}

    legendre(const legendre& o) = default;
    legendre(legendre&& o) = default;

    legendre& operator=(const legendre& o) = default;
    legendre& operator=(legendre&& o) = default;

    void construct_quadrature_rule()
    {
        linalg::symmetric_tridiagonal_matrix<RT, linalg::blas_backend> poly_recurrence(m_w.size(), m_w.size());
        linalg::eigensolver<linalg::symmetric_tridiagonal_matrix<RT, linalg::blas_backend>> eigsolver(m_w.size());
        linalg::vector<RT, linalg::blas_backend> vals(m_w.size());  vals.fill_zeros();
        linalg::matrix<RT, linalg::blas_backend> vecs(m_w.size(), m_w.size());  vecs.fill_zeros();

        for(size_t i=0; i<m_w.size(); ++i)
        {
            size_t k = i+1;
            poly_recurrence(i, i) = 0.0;
            if(i+1 < m_w.size())
            {
                poly_recurrence(i, i+1) = k*sqrt(1.0/((2*k-1)*(2*k+1)));
            }
        }

        eigsolver(poly_recurrence, vals, vecs);
        m_x = vals;
        for(size_t i=0; i<m_w.size(); ++i)
        {
            vals(i) = 2.0*vecs(0, i)*vecs(0, i);
        }
        m_w = vals;
        
    }

    void resize(size_t N)
    {
        if(N != m_w.size())
        {
            m_w.resize(N);
            m_x.resize(N);
            construct_quadrature_rule();
        }
    }

    size_t size() const{return m_w.size();}

    template <typename F>
    auto operator()(F&& f, RT a = -1, RT b = 1) const -> decltype((f(a)))
    {
        RT aa = a;   RT bb = b;
        if(a > b){bb = a;   aa = b;}
        
        return legendre_quad<T, backend>::impl(std::forward<F>(f), aa, bb, m_w, m_x);
    }


    RT weight_function(RT /* x */){return 1.0;}
    
    const linalg::vector<RT, backend>& w() const{return m_w;}
    const linalg::vector<RT, backend>& x() const{return m_x;}

    RT w(size_t i, RT a = -1, RT b = 1) const
    {
        ASSERT(i < m_w.size(), "Index out of bounds.");
        RT aa = a;   RT bb = b;
        if(a > b){bb = a;   aa = b;}

        RT c = (b-a)/2.0;
        return c*m_w(i);
    }

    RT x(size_t i, RT a = -1, RT b = 1) const
    {
        ASSERT(i < m_x.size(), "Index out of bounds.");
        RT aa = a;   RT bb = b;
        if(a > b){bb = a;   aa = b;}

        RT c = (b-a)/2.0;
        RT d = (b+a)/2.0;
        return c*m_x(i)+d;
    }
protected:

    //the quadrature points
    linalg::vector<RT, backend> m_w;
    linalg::vector<RT, backend> m_x;
};

}
}

#endif

