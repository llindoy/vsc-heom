#ifndef HEOM_BATH_CORRELATIONS_HPP
#define HEOM_BATH_CORRELATIONS_HPP

#include <linalg/linalg.hpp>


//an object for handling the generalised heom expansion.  Here we choose the sigma vector such that 
//\sigma_k = \sqrt(|S_k|) = \sqrt(|\sum_l s_{kl}\phi_l(0)).  This ensures that the terms that increase 
//and decrease the number of excitations are the same as those in the Shi renormalisation for exponential 
//bath terms.
template <typename T> 
class correlation_terms
{
public:
    using complex_type = linalg::complex<T>;

public:
    correlation_terms() : m_nmodes(0), m_diagonal_gamma(true){}
    correlation_terms(size_t n) : correlation_terms()       
    {
        CALL_AND_HANDLE(resize(n), "Failed to construct correlation_terms object.");
    }

    template <typename U1, typename U2>
    correlation_terms(U1&& gamma, U2&& s) : correlation_terms()
    {
        CALL_AND_HANDLE(initialise(std::forward<U1>(gamma), std::forward<U2>(s)), "Failed to construct correlation_terms object.");
    }

    template <typename U1, typename U2, typename U3>
    correlation_terms(U1&& gamma, U2&& s, U3&& a) : correlation_terms()
    {
        CALL_AND_HANDLE(initialise(std::forward<U1>(gamma), std::forward<U2>(s), std::forward<U3>(a)), "Failed to construct correlation_terms object.");
    }
    
    correlation_terms(const correlation_terms& o) = default;
    correlation_terms(correlation_terms&& o) = default;

    correlation_terms& operator=(const correlation_terms& o) = default;
    correlation_terms& operator=(correlation_terms&& o) = default;


    template <typename U1>
    void initialise(U1&& gamma)
    {
        CALL_AND_HANDLE(initialise_gamma(std::forward<U1>(gamma)), "Failed to initialise correlation_terms object.");
    }

    template <typename U1, typename U2>
    void initialise(U1&& gamma, U2&& s)
    {
        CALL_AND_HANDLE(initialise_gamma(std::forward<U1>(gamma)), "Failed to initialise correlation_terms object.");
        CALL_AND_HANDLE(set_symmetric(std::forward<U2>(s)), "Failed to initialise correlation_terms object.");
    }

    template <typename U1, typename U2, typename U3>
    void initialise(U1&& gamma, U2&& s, U3&& a)
    {
        CALL_AND_HANDLE(initialise_gamma(std::forward<U1>(gamma)), "Failed to initialise correlation_terms object.");
        CALL_AND_HANDLE(set_symmetric(std::forward<U2>(s)), "Failed to initialise correlation_terms object.");
        CALL_AND_HANDLE(set_asymmetric(std::forward<U3>(a)), "Failed to initialise correlation_terms object.");
    }
 
protected:
    //functions for initialising the gamma object in different ways
    template <typename U>
    void initialise_gamma(const U& gamma)
    {
        resize(1);
        set_gamma(gamma);
    }

    template <typename U>
    void initialise_gamma(const std::vector<U>& gamma)
    {
        size_t n = gamma.size();
        resize(n);

        set_gamma(gamma);
    }

    template <typename U>
    void initialise_gamma(const linalg::vector<U>& gamma)
    {
        size_t n = gamma.size();
        resize(n);
        set_gamma(gamma);
    }

    template <typename U>
    void initialise_gamma(const linalg::matrix<U>& gamma)
    {
        ASSERT(gamma.size(0) == gamma.size(1), "The gamma matrix must be square.");
        size_t n = gamma.size(0);
        resize(n);
        set_gamma(gamma);
    }

    template <typename U>
    void initialise_gamma(const linalg::matrix<U>& gamma, const linalg::matrix<bool>& gamma_con)
    {
        ASSERT(gamma.size(0) == gamma.size(1), "The gamma matrix must be square.");
        ASSERT(gamma_con.size(0) == gamma_con.size(1), "The gamma_con matrix must be square.");
        ASSERT(gamma_con.size(0) == gamma.size(0), "The two matrices must be square.");
        size_t n = gamma.size(0);
        resize(n);
        set_gamma(gamma, gamma_con);
    }

public:
    template <typename U>
    void set_gamma(const U& gamma)
    {
        ASSERT(m_nmodes == 1, "Gamma object is not the correct size.");
        m_diagonal_gamma = true;
        m_gamma_con.fill_value(false);
        m_gammak(0, 0) = gamma;
        m_gamma_con(0, 0) = true;
    }

    //functions for initialising the gamma object in different ways
    template <typename U>
    void set_gamma(const std::vector<U>& gamma)
    {
        ASSERT(gamma.size() == m_nmodes, "Gamma object is not the correct size.");
        m_diagonal_gamma = true;
        m_gamma_con.fill_value(false);
        for(size_t i=0; i < m_nmodes; ++i)
        {
            m_gammak(i, i) = gamma[i];
            m_gamma_con(i, i) = true;
        }
    }

    template <typename U>
    void set_gamma(const linalg::vector<U>& gamma)
    {
        ASSERT(gamma.size() == m_nmodes, "Gamma object is not the correct size.");
        m_diagonal_gamma = true;
        m_gamma_con.fill_value(false);
        for(size_t i=0; i < m_nmodes; ++i)
        {
            m_gammak(i, i) = gamma(i);
            m_gamma_con(i, i) = true;
        }
    }

    template <typename U>
    void set_gamma(const linalg::matrix<U>& gamma)
    {
        ASSERT(gamma.size(0) == gamma.size(1), "The gamma matrix must be square.");
        ASSERT(gamma.size(0) == m_nmodes, "Gamma object is not the correct size.");
        m_diagonal_gamma = false;
        m_gamma_con.fill_value(false);
        for(size_t i=0; i < m_nmodes; ++i)
        {
            for(size_t j=0; j < m_nmodes; ++j)
            {
                m_gammak(i, j) = gamma(i, j);
                m_gamma_con(i, j) = true;
            }
        }
    }

    template <typename U>
    void set_gamma(const linalg::matrix<U>& gamma, const linalg::matrix<bool>& gamma_con)
    {
        ASSERT(gamma.size(0) == m_nmodes, "Gamma object is not the correct size.");
        ASSERT(gamma.size(0) == gamma.size(1), "The gamma matrix must be square.");
        ASSERT(gamma_con.size(0) == gamma_con.size(1), "The gamma_con matrix must be square.");
        ASSERT(gamma_con.size(0) == gamma.size(0), "The two matrices must be square.");
        m_diagonal_gamma = true;
        for(size_t i=0; i < m_nmodes; ++i)
        {
            for(size_t j=0; j < m_nmodes; ++j)
            {
                m_gammak(i, j) = gamma(i, j);
                m_gamma_con(i, j) = gamma_con(i, j);
                if(i != j && gamma_con(i, j)){m_diagonal_gamma = false;}
            }
        }
    }

    template <typename U>
    void set_gamma(size_t i, size_t j, const U& g)
    {
        ASSERT(i < m_nmodes && j < m_nmodes, "Index out of bounds.");
        m_gammak(i, j) = g; 
        m_gamma_con(i, j) = true;
        if(i != j){m_diagonal_gamma = false;}
    }

    void set_gamma_connectivity(const linalg::matrix<bool>& gamma_con)
    {
        ASSERT(gamma_con.size(0) == m_s.size() && gamma_con.size(1) == m_s.size(), "Invalid gamma connectivity array.");
        m_gamma_con = gamma_con;
    }

    void set_gamma_connectivity(size_t i, size_t j, bool gamma_con)
    {
        ASSERT(i == m_s.size() && j == m_s.size(), "Index out of bounds.");
        m_gamma_con(i, j) = gamma_con;
    }

public:
    template <typename U>
    void set_symmetric(const U& s)
    {
        ASSERT(m_nmodes > 0, "cannot set symmetric element if there are no terms.");
        m_s[0] = s;
        m_sp[0] = true;
    }

    template <typename U>
    void set_symmetric(const std::vector<U>& s)
    {
        ASSERT(s.size() == m_nmodes, "s vector is the wrong size.");

        for(size_t i=0; i < m_nmodes; ++i)
        {
            m_s[i] = s[i];
            m_sp[i] = true;
        }
    }

    template <typename U>
    void set_symmetric(const linalg::vector<U>& s)
    {
        ASSERT(s.size() == m_nmodes, "s vector is the wrong size.");

        for(size_t i=0; i < m_nmodes; ++i)
        {
            m_s[i] = s[i];
            m_sp[i] = true;
        }
    }
    
    template <typename U>
    void set_symmetric(const U& s, size_t i)
    {
        ASSERT(i < m_nmodes, "Index out of bounds.");
        m_s(i) = s; 
        m_sp(i) = true;
    }

    template <typename U>
    void set_asymmetric(const U& a)
    {
        ASSERT(m_nmodes > 0, "cannot set asymmetric element if there are no terms.");
       
        m_a[0] = a;
        m_ap[0] = true;
    }
        
    template <typename U>
    void set_asymmetric(const std::vector<U>& a)
    {
        ASSERT(a.size() == m_nmodes, "a vector is the wrong size.");

        for(size_t i=0; i < m_nmodes; ++i)
        {
            m_a[i] = a[i];
            m_ap[i] = true;
        }
    }

    template <typename U>
    void set_asymmetric(const linalg::vector<U>& a)
    {
        ASSERT(a.size() == m_nmodes, "a vector is the wrong size.");

        for(size_t i=0; i < m_nmodes; ++i)
        {
            m_a[i] = a[i];
            m_ap[i] = true;
        }
    }

    template <typename U>
    void set_asymmetric(const U& a, size_t i)
    {
        ASSERT(i < m_nmodes, "Index out of bounds.");
        m_a(i) = a; 
        m_ap(i) = true;
    }

public:
    const complex_type& sk(size_t k) const
    {
        ASSERT(k < m_s.size(), "Failed to access symmetric element. Index out of bounds.");  
        ASSERT(m_sp(k), "Failed to access symmetric element.  Required symmetric element is not present.");
        return m_s(k);
    }

    const complex_type& ak(size_t k) const
    {
        ASSERT(k < m_a.size(), "Failed to access asymmetric element. Index out of bounds.");  
        ASSERT(m_ap(k), "Failed to access asymmetric element.  Required symmetric element is not present.");
        return m_a(k);
    }

    const complex_type& nuk(size_t k) const {ASSERT(k < m_gammak.size(0), "Failed to access element index out of bounds.");  return m_gammak(k, k);}
    const complex_type& gamma(size_t j, size_t k) const
    {
        ASSERT(j < m_gammak.size(0) && k < m_gammak.size(1) && !m_diagonal_gamma,  "Failed to access element index out of bounds.");  return m_gammak(j, k);
    }

    void resize(size_t n)
    {
        m_diagonal_gamma = true;
        m_nmodes = n;
        CALL_AND_HANDLE(m_gammak.resize(n, n), "Failed to excitation number conserving coefficient matrix.");
        CALL_AND_HANDLE(m_gamma_con.resize(n, n), "Failed to excitation number conserving coefficient matrix.");    m_gamma_con.fill_value(false);
 
        CALL_AND_HANDLE(m_s.resize(n), "Failed to resize the symmetric coefficients object.");
        CALL_AND_HANDLE(m_sp.resize(n), "Failed to resize the symmetric coefficients object.");     m_sp.fill_value(false);
        CALL_AND_HANDLE(m_a.resize(n), "Failed to resize the anti symmetric coefficients object.");
        CALL_AND_HANDLE(m_ap.resize(n), "Failed to resize the symmetric coefficients object.");     m_ap.fill_value(false);
    }


public:
    bool nonzero_gamma(size_t i, size_t j) const{ASSERT(i < m_gamma_con.size(0) && j < m_gamma_con.size(1), "Index out of bounds."); return m_gamma_con(i, j);}
    bool has_diagonal_gamma() const {return m_diagonal_gamma;}

    bool has_symmetric_contribution(size_t i) const
    {
        ASSERT(i < m_nmodes, "Index out of bounds.");
        return m_sp[i];
    }
    bool has_asymmetric_contribution(size_t i) const
    {
        ASSERT(i < m_nmodes, "Index out of bounds.");
        return m_ap[i];
    }
    size_t nmodes() const{return m_nmodes;}

protected:
    size_t m_nmodes;
    bool m_diagonal_gamma;
    
    linalg::matrix<complex_type> m_gammak;      //coefficients connecting ados with the same excitation number
    linalg::matrix<bool> m_gamma_con;           //stores whether a given index of the gamma array is occupied

    linalg::vector<complex_type> m_s;           //the symmetric coefficients included in the definition of the \Theta_k object
    linalg::vector<bool> m_sp;                  //stores whether a given index of the symmetric array is present
    linalg::vector<complex_type> m_a;           //the antisymmetric coefficients included in the definition of the \Theta_k object;
    linalg::vector<bool> m_ap;                  //stores whether a given index of the asymmetric array is present
};

#endif

