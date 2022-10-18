#ifndef MODE_DATA_HPP
#define MODE_DATA_HPP

#include <vector>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <utility>
#include <complex>
#include <unordered_set>
#include <tuple>

#include <linalg/linalg.hpp>
#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include <ttns_lib/ttns.hpp>

template <typename T>
class heom_mode_combination
{
protected:
    using real_type = T;
    using complex_type = ttns::complex<T>;
    std::vector<complex_type> m_nuk;
    std::vector<complex_type> m_ck;

    linalg::matrix<size_t> m_nk;
    linalg::matrix<int32_t> m_nkp1;
    linalg::matrix<int32_t> m_nkm1;

    bool m_nu_set;
    bool m_c_set;
    bool m_ind_set;
public:
    heom_mode_combination() :  m_nu_set(false), m_c_set(false), m_ind_set(false) {}
    heom_mode_combination(size_t n) : m_nuk(n), m_ck(n), m_nu_set(false), m_c_set(false), m_ind_set(false) {}

    void resize(size_t n) 
    {
        m_nu_set = false;    m_c_set = false;    m_ind_set = false;
        CALL_AND_HANDLE(m_nuk.resize(n), "Failed to resize harmonic mode combination object.  Failed to resize the frequency array.");
        CALL_AND_HANDLE(m_ck.resize(n), "Failed to resize harmonic mode combination object.  Failed to resize the coupling array.");
        CALL_AND_HANDLE(m_nk.clear(), "Failed to resize harmonic mode combination object.  Failed to clear the topology array.");
        CALL_AND_HANDLE(m_nkp1.clear(), "Failed to resize harmonic mode combination object.  Failed to clear the topology array.");
        CALL_AND_HANDLE(m_nkm1.clear(), "Failed to resize harmonic mode combination object.  Failed to clear the topology array.");
    }

    //functions for setting the inputs necessary
    template <typename R1>
    void set_nuk(const R1& w)
    {
        if(w.size() != m_nuk.size())
        {
            CALL_AND_HANDLE(resize(w.size()), "Failed to set the frequency array of the harmonic mode combination object.");
        }
        for(size_t i=0; i<w.size(); ++i){m_nuk[i] = w[i];}
        m_nu_set = true;
    }

    template <typename R1>
    void set_nuk(R1 w1, R1 w2)
    {
        std::ptrdiff_t nw = w2 - w1;
        if(nw != static_cast<std::ptrdiff_t>(m_nuk.size()))
        {
            CALL_AND_HANDLE(resize(nw), "Failed to set the frequency array of the harmonic mode combination object.");
        }
        
        size_t c = 0;
        for(R1 a = w1; a != w2; ++a, ++c){m_nuk[c] = *a;}
        m_nu_set = true;
    }

    template <typename R1>
    void set_ck(const R1& c)
    {
        if(c.size() != m_ck.size())
        {
            CALL_AND_HANDLE(resize(c.size()), "Failed to set the coupling array of the harmonic mode combination object.");
        }
        for(size_t i=0; i<c.size(); ++i){m_ck[i] = c[i];}
        m_c_set = true;
    }

    template <typename R1>
    void set_ck(R1 c1, R1 c2)
    {
        std::ptrdiff_t nc = c2 - c1;
        if(nc != static_cast<std::ptrdiff_t>(m_ck.size()))
        {
            CALL_AND_HANDLE(resize(nc), "Failed to set the coupling array of the harmonic mode combination object.");
        }
        
        size_t c = 0;
        for(R1 a = c1; a != c2; ++a, ++c){m_ck[c] = *a;}
        m_c_set = true;
    }

    size_t nados(size_t L, size_t K)
    {
        size_t ret = 1;
        for(size_t k = 0; k < K; ++k)
        {
            ret = (ret*(L+k+1))/(k+1);
        }
        return ret;
    }

    //functions for computing the coupling topology required
    size_t construct_basis_topology(size_t L)
    {
        ASSERT(m_nu_set, "Failed to compute the basis topology.   The frequencies of the harmonic modes have not been set.");

        size_t K = m_nuk.size();
        size_t nstates = nados(L, K);

        if(K > 1)
        {
            linalg::vector<size_t> nt(K);
            m_nk.resize(nstates, K);
            populate_state_vector(m_nk, nt, 0, 0, L);
        }
        else
        {
            m_nk.resize(nstates, K);
            for(size_t i=0; i < nstates; ++i)
            {
                m_nk(i, 0) = i;
            }
        }

        //now we need to create a vector of pairs useful for the searching ofvalue arrays
        std::vector<std::vector<size_t>> nk(nstates);   
        for(size_t i=0; i < nstates; ++i)
        {   
            nk[i].resize(K);
            for(size_t j=0; j<K; ++j){nk[i][j] = m_nk(i,j);}
        }

    
        //now we construct the list of lowering operators
        m_nkm1.resize(nstates, K);  m_nkm1.fill_value(-1.0);
        m_nkp1.resize(nstates, K);  m_nkp1.fill_value(-1.0);

        std::vector<size_t> nt(K);
        for(size_t i=0; i < nstates; ++i)
        {
            nt = nk[i];
            for(size_t j=0; j<K; ++j)
            {
                if(nt[j] == 0){m_nkm1(i, j) = -1;}
                else
                {
                    nt[j] = nk[i][j] - 1;
                    auto itr = std::lower_bound(nk.begin(), nk.end(), nt, 
                                    [](const std::vector<size_t>& a, const std::vector<size_t>& b)
                                    {
                                        for(size_t ii=0; ii<a.size(); ++ii)
                                        {
                                            if(a[ii] < b[ii]){return true;}
                                            else if(a[ii] > b[ii]){return false;}
                                        }
                                        return false;
                                    }
                                );
                    if(itr == nk.end()){m_nkm1(i, j) = -1;}
                    else if(*itr != nt){m_nkm1(i, j) = -1;}
                    else{m_nkm1(i, j) = std::distance(nk.begin(), itr);}
                }

                //now we attempt to find the +1
                nt[j] = nk[i][j] + 1;
                auto itr = std::lower_bound(nk.begin(), nk.end(), nt, 
                                [](const std::vector<size_t>& a, const std::vector<size_t>& b)
                                {
                                    for(size_t ii=0; ii<a.size(); ++ii)
                                    {
                                        if(a[ii] < b[ii]){return true;}
                                        else if(a[ii] > b[ii]){return false;}
                                    }
                                    return false;
                                }
                            );
                if(itr == nk.end()){m_nkp1(i, j) = -1;}
                else if(*itr != nt){m_nkp1(i, j) = -1;}
                else{m_nkp1(i, j) = std::distance(nk.begin(), itr);}

                nt[j] = nk[i][j];
            }
        }

        m_ind_set = true;
        return nstates;
    }


    //functions for computing matrices necessary for spin-boson (and polaron transformed spin-boson) dynamics
    template <typename U, typename backend>
    typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type H0(linalg::diagonal_matrix<U, backend>& H)
    {
        ASSERT(m_nu_set && m_ind_set, "Failed to construct H0.  Either the topology or frequency arrays have not been set.");
        size_t nstates = m_nk.shape(0);     size_t N = m_nk.shape(1);
        std::cout << nstates << std::endl;
        H.resize(nstates, nstates);
        linalg::diagonal_matrix<U> _H(nstates, nstates);
        for(size_t i=0; i<nstates; ++i)
        {
            U val = 0.0;
            for(size_t j = 0; j < N; ++j)
            {
                val += m_nuk[j]*static_cast<T>(m_nk(i, j));
            }
            _H(i, i) = val;
        }
        H = _H;
        std::cerr << "H0" << std::endl;
        std::cerr << H << std::endl;
    }

    //functions for computing matrices necessary for spin-boson (and polaron transformed spin-boson) dynamics
    template <typename U, typename backend>
    typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type Hc_commutator(linalg::csr_matrix<U, backend>& H)
    {
        ASSERT(m_nu_set && m_ind_set && m_c_set, "Failed to construct Hc.  Not all of the topology, frequency and coupling arrays have been set.");
        size_t nstates = m_nk.shape(0);     size_t N = m_nk.shape(1);
        linalg::csr_matrix<U> _H;
        _H.resize(nstates, nstates);        H.resize(nstates, nstates);

        //we determine the number of non-zeros in the coupling element, while also setting the rowptr array
        auto rowptr = _H.rowptr();   rowptr[0] = 0;
        size_t nnz = 0;
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t j=0; j<N; ++j)
            {
                nnz += (m_nkm1(i, j) != -1 ? 1 : 0) + (m_nkp1(i, j) != -1 ? 1 : 0);
            }
            rowptr[i+1] = nnz;
        }
        _H.resize(nnz);
        H.resize(nnz);

        //now we can set the colind and buffer values
        size_t counter = 0;
        auto colind = _H.colind();
        auto buffer = _H.buffer();
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t j=0; j<N; ++j)
            {
                if(m_nkm1(i, j) != -1)
                {
                    T n = static_cast<T>(m_nk(i, j));
                    buffer[counter] = sqrt(n/std::abs(m_ck[j]))*linalg::real(m_ck[j]);
                    colind[counter] = m_nkm1(i, j);
                    ++counter;
                }
            }
            
            for(size_t k=0; k<N; ++k)
            {
                size_t j = N-(k+1);
                if(m_nkp1(i, j) != -1)
                {
                    colind[counter] = m_nkp1(i, j);
                    T n = static_cast<T>(m_nk(i, j)+1);
                    buffer[counter] = sqrt(n*std::abs(m_ck[j]));
                    ++counter;
                }
            }
        }
        H = _H;
    }
    
    //functions for computing matrices necessary for spin-boson (and polaron transformed spin-boson) dynamics
    template <typename U, typename backend>
    typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type Hc_anticommutator(linalg::csr_matrix<U, backend>& H)
    {
        ASSERT(m_nu_set && m_ind_set && m_c_set, "Failed to construct Hc.  Not all of the topology, frequency and coupling arrays have been set.");
        size_t nstates = m_nk.shape(0);     size_t N = m_nk.shape(1);
        linalg::csr_matrix<U> _H;
        _H.resize(nstates, nstates);        H.resize(nstates, nstates);  

        //we determine the number of non-zeros in the coupling element, while also setting the rowptr array
        auto rowptr = _H.rowptr();   rowptr[0] = 0;
        size_t nnz = 0;
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t j=0; j<N; ++j)
            {
                nnz += (m_nkp1(i, j) != -1 ? 1 : 0);
            }
            rowptr[i+1] = nnz;
        }
        _H.resize(nnz);
        H.resize(nnz);

        //now we can set the colind and buffer values
        size_t counter = 0;
        auto colind = _H.colind();
        auto buffer = _H.buffer();

        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t k=0; k<N; ++k)
            {
                size_t j = N-(k+1);
                if(m_nkp1(i, j) != -1)
                {
                    colind[counter] = m_nkp1(i, j);
                    T n = static_cast<T>(m_nk(i, j));
                    buffer[counter] = sqrt(n/std::abs(m_ck[j]))*linalg::imag(m_ck[j]);
                    ++counter;
                }
            }
        }
        H = _H;
    }


    template <typename U, typename backend>
    static inline typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type commutator(const linalg::matrix<U, backend>& H, linalg::matrix<U, backend>& L)
    {
        ASSERT(H.shape(0) == H.shape(1), "Can only construct the Liouvillian for a square matrix.");
        size_t N = H.shape(0);
        L.resize(N*N, N*N);
        for(size_t i1 = 0; i1 < N; ++i1)
        {
            for(size_t i2 = 0; i2 < N; ++i2)
            {
                for(size_t j1 = 0; j1 < N; ++j1)
                {
                    for(size_t j2 = 0; j2 < N; ++j2)
                    {
                        L(i1*N+i2, j1*N+j2) = H(i1, j1) * ( i2 == j2  ?  1.0 : 0.0) - (i1 == j1 ? 1.0 : 0.0) * H(j2, i2);
                    }
                }
            }
        }
    }


    template <typename U, typename backend>
    static inline typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type commutator(const linalg::diagonal_matrix<U, backend>& H, linalg::matrix<U, backend>& L)
    {
        ASSERT(H.shape(0) == H.shape(1), "Can only construct the Liouvillian for a square matrix.");
        size_t N = H.shape(0);
        L.resize(N*N, N*N); L.fill_zeros();
        for(size_t i1 = 0; i1 < N; ++i1)
        {
            for(size_t j1 = 0; j1 < N; ++j1)
            {
                L(i1*N+j1, i1*N+j1) = H(i1, i1) - H(j1, j1);
            }
        }
    }

    template <typename U, typename backend>
    static inline typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type anticommutator(const linalg::matrix<U, backend>& H, linalg::matrix<U, backend>& L)
    {
        ASSERT(H.shape(0) == H.shape(1), "Can only construct the Liouvillian for a square matrix.");
        size_t N = H.shape(0);
        L.resize(N*N, N*N);
        for(size_t i1 = 0; i1 < N; ++i1)
        {
            for(size_t i2 = 0; i2 < N; ++i2)
            {
                for(size_t j1 = 0; j1 < N; ++j1)
                {
                    for(size_t j2 = 0; j2 < N; ++j2)
                    {
                        L(i1*N+i2, j1*N+j2) = H(i1, j1) * ( i2 == j2  ?  1.0 : 0.0) + (i1 == j1 ? 1.0 : 0.0) * H(j2, i2);
                    }
                }
            }
        }
    }


    template <typename U, typename backend>
    static inline typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type spre(const linalg::matrix<U, backend>& H, linalg::matrix<U, backend>& L)
    {
        ASSERT(H.shape(0) == H.shape(1), "Can only construct the Liouvillian for a square matrix.");
        size_t N = H.shape(0);
        L.resize(N*N, N*N);
        for(size_t i1 = 0; i1 < N; ++i1)
        {
            for(size_t i2 = 0; i2 < N; ++i2)
            {
                for(size_t j1 = 0; j1 < N; ++j1)
                {
                    for(size_t j2 = 0; j2 < N; ++j2)
                    {
                        L(i1*N+i2, j1*N+j2) = H(i1, j1) * ( i2 == j2  ?  1.0 : 0.0);
                    }
                }
            }
        }
    }

    template <typename U, typename backend>
    static inline typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type spost(const linalg::matrix<U, backend>& H, linalg::matrix<U, backend>& L)
    {
        ASSERT(H.shape(0) == H.shape(1), "Can only construct the Liouvillian for a square matrix.");
        size_t N = H.shape(0);
        L.resize(N*N, N*N);
        for(size_t i1 = 0; i1 < N; ++i1)
        {
            for(size_t i2 = 0; i2 < N; ++i2)
            {
                for(size_t j1 = 0; j1 < N; ++j1)
                {
                    for(size_t j2 = 0; j2 < N; ++j2)
                    {
                        L(i1*N+i2, j1*N+j2) = (i1 == j1 ? 1.0 : 0.0) * H(j2, i2);
                    }
                }
            }
        }
    }


    template <typename U, typename backend>
    static inline typename std::enable_if<std::is_same<T, U>::value || std::is_same<complex_type, U>::value, void>::type kron(const linalg::matrix<U, backend>& H, const linalg::matrix<U, backend>& H2,linalg::matrix<U, backend>& L)
    {
        ASSERT(H.shape(0) == H.shape(1), "Can only construct the Liouvillian for a square matrix.");
        size_t N = H.shape(0);
        L.resize(N*N, N*N);
        for(size_t i1 = 0; i1 < N; ++i1)
        {
            for(size_t i2 = 0; i2 < N; ++i2)
            {
                for(size_t j1 = 0; j1 < N; ++j1)
                {
                    for(size_t j2 = 0; j2 < N; ++j2)
                    {
                        L(i1*N+i2, j1*N+j2) = H(i1, j1) * H2(j2, i2);
                    }
                }
            }
        }
    }


    template <typename backend>
    static inline void onsite(const linalg::diagonal_matrix<complex_type, backend>& H, const linalg::matrix<complex_type, backend>& F, const T& delta, linalg::matrix<complex_type, backend>& L)
    {
        commutator(H, L);       

        linalg::matrix<complex_type, backend> Fcomm2;
        {
            linalg::matrix<complex_type, backend> Fcomm;
            commutator(F, Fcomm); 
            linalg::matrix<complex_type, backend> Fcommt(Fcomm);
            Fcomm2  = Fcommt*Fcomm;
        }

        //L += complex_type(0, 1)*delta*Fcomm2;
    }
      
    //the indexing here is completely wrong
    template <typename backend>
    void heom_operator(const linalg::diagonal_matrix<complex_type, backend>& H, const linalg::matrix<complex_type, backend>& F, const T& delta, linalg::csr_matrix<complex_type, backend>& A)
    {
        std::cerr << "coeffs" << std::endl;
        for(size_t i = 0; i < m_ck.size(); ++i)
        {
            std::cerr << sqrt(m_ck[i]) << " " << m_nuk[i] << std::endl;
        }
        size_t nhilb = H.shape(0);
        linalg::matrix<complex_type, backend> L;
        commutator(H, L);       
        linalg::matrix<complex_type, backend> Q = F;


        {
            linalg::matrix<complex_type, backend> Qd = linalg::adjoint(F);
            linalg::matrix<complex_type, backend> F2 = Qd*Q;

            std::cerr << Q << " " << Qd << std::endl;

            linalg::matrix<complex_type, backend> op;
            linalg::matrix<complex_type, backend> Fop;
            linalg::matrix<complex_type, backend> Fop2;

            spre(F2, Fop);      op = Fop;
            spost(F2, Fop);     op += Fop;

            spre(Q, Fop);      
            spost(Qd, Fop2);    op -= 2.0*Fop*Fop2;

            L = 100.0*L + complex_type(0, 1)*delta*op;
        }


        size_t nliouv = L.shape(0);

        ASSERT(m_nu_set && m_ind_set && m_c_set, "Failed to construct Hc.  Not all of the topology, frequency and coupling arrays have been set." );
        size_t nstates = m_nk.shape(0);     size_t N = m_nk.shape(1);
        linalg::csr_matrix<complex_type> _H;
        _H.resize(nstates*nliouv, nstates*nliouv);        A.resize(nstates*nliouv, nstates*nliouv);

        //we determine the number of non-zeros in the coupling element, while also setting the rowptr array
        auto rowptr = _H.rowptr();   rowptr[0] = 0;
        size_t nnz = 0;
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t l=0; l<nliouv; ++l)
            {
                for(size_t j=0; j<N; ++j)
                {
                    nnz += ((m_nkm1(i, j) != -1 ? 1 : 0) + (m_nkp1(i, j) != -1 ? 1 : 0))*(2*nhilb-1);
                }
                nnz += nliouv;
                rowptr[i*nliouv+l+1] = nnz;
            }
        }
        _H.resize(nnz);
        A.resize(nnz);

        //now we can set the colind and buffer values
        size_t counter = 0;
        auto colind = _H.colind();
        auto buffer = _H.buffer();
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t hi1 = 0; hi1 < nhilb; ++hi1)
            {
                for(size_t hi2 = 0; hi2 < nhilb; ++hi2)
                {
                    //add on the lower diagonal terms
                    for(size_t j=0; j<N; ++j)
                    {
                        if(m_nkm1(i, j) != -1)
                        {
                            for(size_t hj1 = 0; hj1 < nhilb; ++hj1)
                            {
                                for(size_t hj2 = 0; hj2 < nhilb; ++hj2)
                                {
                                    if(hi1 == hj1 || hi2 == hj2)
                                    {
                                        T n = static_cast<T>(m_nk(i, j));
                                        colind[counter] = m_nkm1(i, j)*nliouv+(hj1*nhilb + hj2);
                                        buffer[counter] = -sqrt(n/std::abs(m_ck[j]))*
                                                                (m_ck[j]*Q(hi1, hj1)*(hi2 == hj2 ? 1.0 : 0.0) - linalg::conj(m_ck[j])*(hi1 == hj1 ? 1.0 : 0.0)*Q(hj2, hi2));
                                        //buffer[counter] = n*(m_ck[j]*Q(hi1, hj1)*(hi2 == hj2 ? 1.0 : 0.0) - linalg::conj(m_ck[j])*(hi1 == hj1 ? 1.0 : 0.0)*Q(hj2, hi2));
                                        ++counter;
                                    }
                                }
                            }
                        }
                    }
                    //add on the diagonal in ado space terms
                    size_t li = hi1*nhilb+hi2;
                    for(size_t lj = 0; lj < nliouv; ++lj)
                    {
                        buffer[counter] = L(li, lj);
                        colind[counter] = i*nliouv+lj;
                        if(li == lj)
                        {
                            complex_type val = 0.0;
                            for(size_t j = 0; j < N; ++j)
                            {
                                val += m_nuk[j]*static_cast<T>(m_nk(i, j));
                            }
                            buffer[counter] += complex_type(0, -1)*val;
                        }
                        ++counter;
                    }
                    //add on the upper diagonal terms
                    for(size_t k=0; k<N; ++k)
                    {
                        size_t j = N-(k+1);
                        if(m_nkp1(i, j) != -1)
                        {
                            for(size_t hj1 = 0; hj1 < nhilb; ++hj1)
                            {
                                for(size_t hj2 = 0; hj2 < nhilb; ++hj2)
                                {
                                    if(hi1 == hj1 || hi2 == hj2)
                                    {

                                        T n = static_cast<T>(m_nk(i, j)+1.0);
                                        colind[counter] = m_nkp1(i, j)*nliouv+(hj1*nhilb + hj2);
                                        buffer[counter] = -sqrt(n*std::abs(m_ck[j]))*(Q(hi1, hj1)*(hi2 == hj2 ? 1.0 : 0.0) - (hi1 == hj1 ? 1.0 : 0.0)*Q(hj2, hi2));
                                        //buffer[counter] = (Q(hi1, hj1)*(hi2 == hj2 ? 1.0 : 0.0) - (hi1 == hj1 ? 1.0 : 0.0)*Q(hj2, hi2));
                                        ++counter;
                                    }
                                }
                            }
                        }
                    }
                    std::cerr << i*nliouv+hi1*nhilb+hi2+1 << " " << rowptr[i*nliouv+hi1*nhilb+hi2+1] << " " << counter << std::endl;
                    //ASSERT(rowptr[i*nliouv+hi1*nhilb+hi2+1] == counter, "Something went wrong.");
                }
            }
        }
        A = _H;
    }

    template <typename backend>
    void heom_coupling_operator(const linalg::matrix<complex_type, backend>& F, linalg::csr_matrix<complex_type, backend>& Fcomm/*,  linalg::csr_matrix<complex_type, backend>& Fanticomm  */)
    {
        size_t nhilb = F.shape(0);
        size_t nliouv = nhilb*nhilb;

        ASSERT(m_ind_set, "Failed to construct Hc.  Not all of the topology, frequency and coupling arrays have been set.");
        size_t nstates = m_nk.shape(0);     size_t N = m_nk.shape(1);
        linalg::csr_matrix<complex_type> _H;
        _H.resize(nstates*nliouv, nstates*nliouv);        Fcomm.resize(nstates*nliouv, nstates*nliouv);     //Fanticomm.resize(nstates*nliouv, nstates*nliouv);

        //we determine the number of non-zeros in the coupling element, while also setting the rowptr array
        auto rowptr = _H.rowptr();   rowptr[0] = 0;
        size_t nnz = 0;
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t l=0; l<nliouv; ++l)
            {
                nnz += (2*nhilb-1);
                rowptr[i*nliouv+l+1] = nnz;
            }
        }
        _H.resize(nnz);
        Fcomm.resize(nnz);

        //now we can set the colind and buffer values
        size_t counter = 0;
        auto colind = _H.colind();
        auto buffer = _H.buffer();
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t hi1 = 0; hi1 < nhilb; ++hi1)
            {
                for(size_t hi2 = 0; hi2 < nhilb; ++hi2)
                {
                    for(size_t hj1 = 0; hj1 < nhilb; ++hj1)
                    {
                        for(size_t hj2 = 0; hj2 < nhilb; ++hj2)
                        {
                            size_t lj = hj1*nhilb+hj2;
                            if(hi1 == hj1 || hi2 == hj2)
                            {
                                buffer[counter] = (F(hi1, hj1)*(hi2 == hj2 ? 1.0 : 0.0) - (hi1 == hj1 ? 1.0 : 0.0)*F(hj2, hi2));
                                colind[counter] = i*nliouv+lj;
                                ++counter;
                            }
                        }
                    }
                }
            }
        }
        Fcomm = _H;

        /*  
        counter = 0;
        for(size_t i=0; i<nstates; ++i)
        {
            for(size_t hi1 = 0; hi1 < nhilb; ++hi1)
            {
                for(size_t hi2 = 0; hi2 < nhilb; ++hi2)
                {
                    for(size_t hj1 = 0; hj1 < nhilb; ++hj1)
                    {
                        for(size_t hj2 = 0; hj2 < nhilb; ++hj2)
                        {
                            size_t lj = hj1*nhilb+hj2;
                            if(hi1 == hj1 || hi2 == hj1)
                            {
                                buffer[counter] = (F(hi1, hj1)*(hi2 == hj2 ? 1.0 : 0.0) + (hi1 == hj1 ? 1.0 : 0.0)*F(hj2, hi2));
                                colind[counter] = i*nliouv+lj;
                                ++counter;
                            }
                        }
                    }
                }
            }
        }
        Fanticomm = _H;
        */
    }
protected:
    size_t populate_state_vector(linalg::matrix<size_t>& nk, linalg::vector<size_t>& nt, size_t curr_ind, size_t mode_ind, size_t L)
    {
        if(mode_ind < m_nuk.size())
        {
            for(size_t i=0; i<=L; ++i)
            {   
                nt[mode_ind] = i;
                curr_ind = populate_state_vector(nk, nt, curr_ind, mode_ind+1, L-i);
            }
        }
        else
        {   
            std::cerr << curr_ind << ": " << std::endl;
            for(size_t i=0; i < m_nuk.size(); ++i){std::cerr << nt[i] << " ";} std::cerr << std::endl;
            nk[curr_ind] = nt;  ++curr_ind;
        }
        return curr_ind;
    }
};


#endif

