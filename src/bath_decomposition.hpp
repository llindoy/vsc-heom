#ifndef HEOM_BATH_DECOMPOSITION_HPP
#define HEOM_BATH_DECOMPOSITION_HPP

template <typename T> 
class bose_decomposition
{
public:
    bose_decomposition() : m_K(0) {}
    bose_decomposition(size_t K, T beta) : m_K(K), m_beta(beta), m_initialised(false) {}
    virtual ~bose_decomposition(){}

    const T& beta() const{return m_beta;}
    T& beta(){return m_beta;}

    const size_t& K() const{return m_K;}
    size_t& K(){return m_K;}

    virtual T nu(size_t i) = 0;
    virtual T eta(size_t i) = 0;
protected:
    bool m_initialised;
    size_t m_K;
    T m_beta;
};


template <typename T> 
class matsubara_decomposition : public bose_decomposition<T>
{
public:
    matsubara_decomposition() : bose_decomposition<T>() {}
    matsubara_decomposition(size_t K, T beta) : bose_decomposition<T>(K, beta) {}

    ~matsubara_decomposition(){}

    T nu(size_t i)
    {
        return 2*M_PI*(i+1)/(this->m_beta);
    }
    T eta(size_t i) 
    {
        return 1.0;
    }
};


template <typename T> 
class pade_decomposition : public bose_decomposition<T>
{
public:
    pade_decomposition() : bose_decomposition<T>() {}
    pade_decomposition(size_t K, T beta) : bose_decomposition<T>(K, beta) 
    {
        init();
    }

    ~pade_decomposition(){}

    T nu(size_t i)
    {
        if(this->m_initialised)
        {
            init();
        }
        ASSERT(i < m_eta.size(), "Index out of bounds.");
        return m_nu[i];
    }
    T eta(size_t i)
    {
        if(this->m_initialised)
        {
            init();
        }
        ASSERT(i < m_eta.size(), "Index out of bounds.");
        return m_eta[i];
    }


    void init()
    {
        size_t N = this->m_K;
        T beta = this->m_beta;
        m_eta.resize(N);
        m_nu.resize(N);
        linalg::vector<T> v
        (
            2*N, 
            [](size_t i)
            {
                T m = i+1.0;
                return 1.0/(std::sqrt((2*m+1.0)*(2*(m+1.0)+1.0)));
            }
        );
        linalg::symmetric_tridiagonal_matrix<T> L(2*N, 2*N);  
        
        for(size_t i = 0; i < 2*N; ++i){L(i, i)= 0;}
        for(size_t i=0; i < 2*N-1; ++i) 
        {
            L(i, i+1) = v(i);
            L(i+1, i) = v(i);
        }

        linalg::eigensolver<linalg::symmetric_tridiagonal_matrix<T>> solver(2*N);
        linalg::diagonal_matrix<T> c1;
        linalg::vector<T> xi(N);
        solver(L, c1);
        for(size_t i = N; i < c1.size(); ++i){xi(i-N) = 2.0/c1(i, i);}

        v.fill
        (
            [](size_t i)
            {
                T m = i+1.0;
                return 1.0/(std::sqrt((2*m+3.0)*(2*(m+1.0)+3.0)));
            }
        );
        linalg::symmetric_tridiagonal_matrix<T> Lp(2*N-1, 2*N-1);  

        for(size_t i = 0; i < 2*N-1; ++i){Lp(i, i)= 0;}
        for(size_t i=0; i < 2*N-2; ++i) 
        {
            Lp(i, i+1) = v(i);
            Lp(i+1, i) = v(i);
        }

        linalg::diagonal_matrix<T> c2;
        linalg::vector<T> zeta(N-1);
        solver(Lp, c2);

        for(size_t i = N; i < c2.size(); ++i){zeta(i-N) = 2.0/c2(i, i);}

        for(size_t j = 0; j < N; ++j)   
        {
            T prod = 1;
            for(size_t k=0; k< N-1; ++k)
            {
                if(k == j)
                {
                    prod *= (zeta[k]*zeta[k] -xi[j]*xi[j]) / (xi[N-1]*xi[N-1] - xi[j]*xi[j]);
                }
                else    
                {
                    prod *= (zeta[k]*zeta[k] - xi[j]*xi[j]) / (xi[k]*xi[k]-xi[j]*xi[j]);
                }
            }

            m_eta[j] = (N*N+1.5*N)*prod;
            m_nu[j] = xi[j]/beta;
        }
    }
protected:
    linalg::vector<T> m_nu;
    linalg::vector<T> m_eta;
};

#endif

