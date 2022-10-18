#ifndef HEOM_OPERATOR_LOW_MEM_HPP
#define HEOM_OPERATOR_LOW_MEM_HPP

#include <linalg/linalg.hpp>

#include "occupation_number_basis_indexing.hpp"
#include "heom_bath_terms.hpp"
#include "superoperator_utilities.hpp"


//an object for constructing generic HEOM bath operators.  For an arbitrary spectral density the generalised heom can be written in a system bath form 
//with an operator acting on the bath degrees of freedom (no action on the system degrees of freedom) and up to two terms that couple the system and bath degrees
//of freedom.  A term accounting for the symmetric part of the bath correlation function and potentially a term accounting for the asymmetric part.  
template <typename T>
class heom_bath_operator_low_mem
{
public:
    using complex_type = linalg::complex<T>;

protected:
    truncated_occupation_number_basis m_indexing;

    size_t m_K;
    size_t m_nhilb;
    bool m_diagonal_bath;
    bool m_has_asymmetric_coupling;
    bool m_interactions_initialised;

    bool m_has_truncation_contribution;
    bool m_has_non_constant_truncation_contribution;
    bool m_has_constant_truncation_contribution;

    bool m_indexing_initialised;
    bool m_rescale_couplings;
    std::vector<complex_type> m_nuk;
    std::vector<size_t> m_jn;

    linalg::matrix<complex_type> m_Ls;
        
    linalg::vector<linalg::diagonal_matrix<complex_type>> m_Lkp;
    linalg::vector<linalg::diagonal_matrix<complex_type>> m_Lkm;

public: 
    heom_bath_operator_low_mem() : m_K(0), m_diagonal_bath(false), m_has_asymmetric_coupling(false), m_interactions_initialised(false), m_indexing_initialised(false), m_rescale_couplings(true) {}

    heom_bath_operator_low_mem(const heom_bath<T>& o, const size_t& L) : heom_bath_operator_low_mem()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }
    heom_bath_operator_low_mem(const std::vector<heom_bath<T>>& o, const size_t& L) : heom_bath_operator_low_mem()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }
    heom_bath_operator_low_mem(const std::vector<std::vector<heom_bath<T>>>& o, const size_t& L) : heom_bath_operator_low_mem()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }

    heom_bath_operator_low_mem(const heom_bath_operator_low_mem& o) = default;
    heom_bath_operator_low_mem(heom_bath_operator_low_mem&& o) = default;

    heom_bath_operator_low_mem& operator=(const heom_bath_operator_low_mem& o) = default;
    heom_bath_operator_low_mem& operator=(heom_bath_operator_low_mem&& o) = default;

    const bool& rescale_coupling_constants() const{return m_rescale_couplings;}
    bool& rescale_coupling_constants(){return m_rescale_couplings;}

    void initialise(const size_t& L)
    {
        CALL_AND_RETHROW(initialise_indexing(L));
    }

    void initialise_cutoff(const T& wmax)
    {
        CALL_AND_RETHROW(initialise_indexing(wmax));
    }

    void initialise_cutoff(const T& wmax, const T& maxfactor)
    {
        CALL_AND_RETHROW(initialise_indexing(wmax, maxfactor));
    }

    void initialise_indexing(const size_t& L)
    {
        m_indexing.clear();
        std::cerr << L << " " << m_K << std::endl;
        CALL_AND_HANDLE(m_indexing.initialise(L, m_K), "Failed to initialise heom operator object.  Failed to initialise the basis indexing object.");
        std::cerr << L << " " << m_K << std::endl;
        m_indexing_initialised = true;
    }

    void initialise_indexing(const T& wmax)
    {
        std::vector<T> nuk(m_K);
        for(size_t i=0; i < m_K; ++i){nuk[i] = linalg::real(m_nuk[i]);}
        
    
        m_indexing.clear();
        CALL_AND_HANDLE(m_indexing.initialise(nuk, wmax), "Failed to initialise heom operator object.  Failed to initialise the basis indexing object.");
        m_indexing_initialised = true;
    }

    void initialise_indexing(const T& wmax, const T& maxfactor)
    {
        std::vector<T> nuk(m_K);
        for(size_t i=0; i < m_K; ++i)
        {
            nuk[i] = linalg::real(m_nuk[i]);
            if(nuk[i] > maxfactor*nuk[0]){nuk[i] = maxfactor*nuk[0];}
        }
    
        m_indexing.clear();
        CALL_AND_HANDLE(m_indexing.initialise(nuk, wmax), "Failed to initialise heom operator object.  Failed to initialise the basis indexing object.");
        m_indexing_initialised = true;
    }

    void clear()
    {
        m_indexing.clear();
        m_K = 0;
        m_diagonal_bath = true;
        m_has_asymmetric_coupling = false;
        m_has_truncation_contribution = false;
        m_has_non_constant_truncation_contribution = false;
        m_has_constant_truncation_contribution = false;
        m_indexing_initialised = false;
        m_rescale_couplings = false;
        m_nuk.clear();

        m_Ls.clear();
        m_Lkp.clear();
        m_Lkm.clear();
        m_jn.clear();
    }



    void set_interactions(std::shared_ptr<superoperator::superoperator<complex_type>> H, const linalg::diagonal_matrix<complex_type>& S, const heom_bath<T>& o)
    {
        try
        {
            std::vector<linalg::diagonal_matrix<complex_type>> _S(1);
            _S[0] = S;

            std::vector<heom_bath<T>> _o(1);    _o[0] = o;

            CALL_AND_RETHROW(set_interactions(H, _S, _o));
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to initialise sparse heom object.");
        }
    }

    void set_interactions(std::shared_ptr<superoperator::superoperator<complex_type>> H, const std::vector<linalg::diagonal_matrix<complex_type>>& S, const std::vector<heom_bath<T>>& o)
    {
        ASSERT(S.size() == o.size(), "The number of bath correlation terms and system contributions to the system-bath coupling operators is not consistent.");
        size_t K = 0;

        //first determine the total number of modes
        for(size_t i=0; i < o.size(); ++i)
        {
            for(size_t ti = 0; ti < o[i].size(); ++ti)
            {
                K += o[i][ti].nmodes();
                ASSERT(o[i][ti].has_diagonal_gamma(), "sparse heom currently only support standard HEOM operators.");
            }
        }

        //now we set up the system operator
        CALL_AND_HANDLE(system_operator(H), "Failed to construct system operator.");

        //now set up the frequency array and the two coupling operator terms
        m_nuk.resize(K);
        m_Lkp.resize(K);
        m_Lkm.resize(K);
        m_jn.resize(K);
        size_t counter = 0;
        m_K = K;
        m_nhilb = H->hilb_dim();
        linalg::diagonal_matrix<complex_type> mS;
        for(size_t i=0; i < o.size(); ++i)
        {
            for(size_t ti = 0; ti < o[i].size(); ++ti)
            {
                for(size_t j=0; j < o[i][ti].nmodes(); ++j)
                {
                    m_nuk[counter] = o[i][ti].nuk(j);

                    CALL_AND_HANDLE(form_Lkp(S[i], o[i], ti, j, mS, m_Lkp[counter]), "Failed to construct L_k+ operator.");
                    CALL_AND_HANDLE(form_Lkm(S[i], o[i], ti, j, mS, m_Lkm[counter]), "Failed to construct L_k- operator.");
                    
                    ASSERT(m_Lkp[counter].shape(0) == m_Ls.shape(0) && m_Lkp[counter].shape(1) == m_Ls.shape(1), "The Lkp operator formed does not have the same size as the Ls operator.");
                    ASSERT(m_Lkm[counter].shape(0) == m_Ls.shape(0) && m_Lkm[counter].shape(1) == m_Ls.shape(1), "The Lkm operator formed does not have the same size as the Ls operator.");
        
                    ++counter;
                }
            }
        }

        m_indexing_initialised = false;
        m_interactions_initialised = true;
    }

    size_t nados() const{return m_indexing.nstates();}
    size_t nmodes() const{return m_K;}

public:
    /*
     * Function for constructing the general HEOM evolution operator for a system interacting with a set of
     * baths with differing system coupling operators.  Here the user must specify as many system operators
     * as there are distinct terms stored by this object.
     */

    void apply(const linalg::vector<complex_type>& A, linalg::vector<complex_type>& B, complex_type coeff = complex_type(1, 0))
    {
        try
        {
            ASSERT(m_interactions_initialised && m_indexing_initialised, "Need to setup object before using.");

            size_t nhilb = m_nhilb;
            size_t _nados = nados();
            m_jn.resize(m_K);

            B.fill_zeros();
            auto ados = A.reinterpret_shape(_nados, nhilb*nhilb);
            auto Bados = B.reinterpret_shape(_nados, nhilb*nhilb);

            for(size_t iado=0; iado < _nados; ++iado)
            {
                auto Bi = Bados[iado];
                m_indexing.state_at(iado, m_jn);
                complex_type gamma_n = complex_type(0, -1)*diagonal_bath_term(m_jn);

                //apply t
                for(size_t k=0; k < m_K; ++k)
                {
                    if(m_indexing._contains_lowered_state(iado, k))
                    {
                        size_t ind = m_indexing._get_lowered_index(iado, k);
                        Bi += (sqrt(m_jn[k]*1.0))*m_Lkp[k]*ados[ind];
                    }
                }
                
                //apply the diagonal term
                {
                    Bi += gamma_n * ados[iado];
                    Bi += m_Ls*ados[iado];
                }


                //now we add the terms that couple this ado to those with larger numbers of excitations
                for(size_t k=0; k < m_K; ++k)
                {
                    if(m_indexing._contains_raised_state(iado, k))
                    {
                        size_t ind = m_indexing._get_raised_index(iado, k);
                        Bi += sqrt(m_jn[k]+1.0)*m_Lkm[k]*ados[ind];
                    }
                }
            }
            B*=coeff;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to apply HEOM evolution operator.");
        }
    } 

protected:
    complex_type diagonal_bath_term(const std::vector<size_t>& n) const
    {
        complex_type ci(0, 0);
        for(size_t k = 0; k < m_K; ++k){ci += m_nuk[k] * static_cast<T>(n[k]);}
        return ci;
    }

protected:
    void form_Lkp(const linalg::diagonal_matrix<complex_type>& S, const heom_bath<T>& o, size_t ti, size_t k, linalg::diagonal_matrix<complex_type>& Ms, linalg::diagonal_matrix<complex_type>& M)
    {
        try
        {
            size_t nhilb = S.shape(0);
            M.resize(nhilb*nhilb, nhilb*nhilb);
            //and scale it by the appropriate scaling factor
            T prefactor = 1.0;

            prefactor = std::sqrt(linalg::abs(linalg::real(o[ti].sk(k))));
            for(size_t i =0; i  < nhilb; ++i)
            {
                for(size_t j = 0; j < nhilb; ++j)
                {
                    if(o[ti].has_symmetric_contribution(k))
                    {
                        M(i*nhilb+j,i*nhilb+j) = (o[ti].sk(k)/prefactor)*(S(i,i) - S(j,j));
                    }
                    //now we add 
                    if(o[ti].has_asymmetric_contribution(k))
                    {                
                        M(i*nhilb+j, i*nhilb+j) += ( complex_type(0, 1)*o[ti].ak(k)/prefactor)*(S(i, i)+S(j,j));
                    }
                }
            }
            
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct Lk+ operator.");
        }
    }

    void form_Lkm(const linalg::diagonal_matrix<complex_type>& S, const heom_bath<T>& o, size_t ti, size_t k, linalg::diagonal_matrix<complex_type>& Ms, linalg::diagonal_matrix<complex_type>& M)
    {
        try
        {
            size_t nhilb = S.shape(0);
            M.resize(nhilb*nhilb, nhilb*nhilb);
            //and scale it by the appropriate scaling factor
            T prefactor = 1.0;
            //and scale it by the appropriate scaling factor
            prefactor = std::sqrt(linalg::abs(linalg::real(o[ti].sk(k))));

            for(size_t i =0; i  < nhilb; ++i)
            {
                for(size_t j = 0; j < nhilb; ++j)
                {
                    M(i*nhilb+j,i*nhilb+j) = prefactor*(S(i,i) - S(j,j));
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct Lk- operator.");
        }
    }
    

    void system_operator(std::shared_ptr<superoperator::superoperator<complex_type>> H)
    {
        try
        {
            linalg::csr_matrix<complex_type> L;
            H->commutator(L);

            superoperator::csr_to_dense(L, m_Ls);

            //linalg::matrix<complex_type> So;
            //for(size_t tt = 0; tt < S.size(); ++tt)
            //{
            //    if(o[tt].has_constant_truncation())
            //    {
            //        o[tt].truncation()->apply(S[tt], m_Ls);
            //    }
            //}
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct system operator matrix.");
        }
    }
};

#endif

