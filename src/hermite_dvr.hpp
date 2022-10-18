#ifndef HERMITE_DVR_HPP
#define HERMITE_DVR_HPP

#include <linalg/linalg.hpp>
#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "utils/orthopol.hpp"

//sets up the hermite dvr
template <typename T>
class hermite_dvr
{
    using complex_type = T;
public:
    hermite_dvr() : m_N(0),m_w0(1),  m_quadrature_constructed(false), m_polynomials_constructed(false){}
    hermite_dvr(size_t N) : m_N(N),m_w0(1),  m_quadrature_constructed(false), m_polynomials_constructed(false){}
    hermite_dvr(size_t N, const T& w0) : m_N(N),m_w0(w0), m_quadrature_constructed(false), m_polynomials_constructed(false){}

    const size_t& Ndvr() const{return m_N;}
    size_t& Ndvr()
    {
        m_polynomials_constructed = false;
        m_quadrature_constructed = false;
        return m_N;
    }

    const T& w0() const{return m_w0;}
    T& w0(){return m_w0;}

    const linalg::vector<T>& x() const{ASSERT(m_quadrature_constructed && m_polynomials_constructed, "Nodes have not been constructed."); return m_poly.nodes();}
    const linalg::vector<T>& w() const{ASSERT(m_quadrature_constructed && m_polynomials_constructed, "Nodes have not been constructed."); return m_poly.weights();}

    void d2dq2(linalg::matrix<complex_type>& k)
    {
        if(!m_polynomials_constructed){construct_polynomials();}
        if(!m_quadrature_constructed){construct_quadrature();}
        size_t n = m_N;
        for(size_t i = 0; i < n; ++i)
        {
            T xi = m_poly.nodes()(i);
            for(size_t j = 0; j < n; ++j)
            {
                T xj = m_poly.nodes()(j);
                if(i == j){k(i, j) = -2.0*(n-1)/3.0 - 0.5 + xi*xi/3.0;}
                else
                {
                    k(i, j) = ((i+j) % 2 == 0 ? 1.0 : -1.0)*(0.5 - 2.0/((xi-xj)*(xi-xj)));
                }
            }
        }
    }

    T V(const T& x)
    {
        return x*x;
    }

    void H0(linalg::matrix<complex_type>& H)
    {
        if(!m_polynomials_constructed){construct_polynomials();}
        if(!m_quadrature_constructed){construct_quadrature();}
        size_t n = m_N;

        d2dq2(H);
        H *= -1.0;
        for(size_t i = 0; i < n; ++i)
        {
            H(i, i) += V(m_poly.nodes()(i));
        }
        H *= m_w0/2.0;
        for(size_t i = 0; i < n; ++i)
        {
            H(i, i) -= m_w0/2.0;
        }
    }

    void construct_quadrature()
    {
        if(!m_polynomials_constructed)
        {
            construct_polynomials();
        }
        if(!m_quadrature_constructed)
        {
            m_poly.compute_nodes_and_weights(); 
            m_quadrature_constructed = true;
        }
    }
    void construct_polynomials()
    {
        CALL_AND_HANDLE(eos::hermite_polynomial(m_poly, m_N), "Failed to constructe Hermite dvr object.");
        m_quadrature_constructed = false;
        m_polynomials_constructed = true;
    }

    const eos::orthopol<T>& poly() const{return m_poly;}
protected:
    eos::orthopol<T> m_poly;
    
    size_t m_N;
    T m_w0;
    bool m_quadrature_constructed;
    bool m_polynomials_constructed;
};

#endif

