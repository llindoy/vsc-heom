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
class heom_bath_operator_sparse
{
public:
    using complex_type = linalg::complex<T>;

protected:
    truncated_occupation_number_basis m_indexing;

    std::vector<size_t> m_Kb;
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
    linalg::vector<complex_type> m_c1k;
    linalg::vector<complex_type> m_c2k;
    linalg::vector<complex_type> m_c3k;
    std::vector<size_t> m_jn;

    linalg::matrix<complex_type> m_Hs;
    std::vector<linalg::matrix<complex_type>> m_S;
    std::vector<heom_bath<T>> m_terms;

    linalg::matrix<complex_type> m_Srho;
    linalg::matrix<complex_type> m_rhoS;

public: 
    heom_bath_operator_sparse() : m_K(0), m_diagonal_bath(false), m_has_asymmetric_coupling(false), m_interactions_initialised(false), m_indexing_initialised(false), m_rescale_couplings(true) {}

    heom_bath_operator_sparse(const heom_bath<T>& o, const size_t& L) : heom_bath_operator_sparse()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }
    heom_bath_operator_sparse(const std::vector<heom_bath<T>>& o, const size_t& L) : heom_bath_operator_sparse()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }
    heom_bath_operator_sparse(const std::vector<std::vector<heom_bath<T>>>& o, const size_t& L) : heom_bath_operator_sparse()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }

    heom_bath_operator_sparse(const heom_bath_operator_sparse& o) = default;
    heom_bath_operator_sparse(heom_bath_operator_sparse&& o) = default;

    heom_bath_operator_sparse& operator=(const heom_bath_operator_sparse& o) = default;
    heom_bath_operator_sparse& operator=(heom_bath_operator_sparse&& o) = default;

    const bool& rescale_coupling_constants() const{return m_rescale_couplings;}
    bool& rescale_coupling_constants(){return m_rescale_couplings;}

    void initialise(const size_t& L)
    {
        CALL_AND_RETHROW(initialise_common());
        CALL_AND_RETHROW(initialise_indexing(L));
    }

    void initialise_cutoff(const T& wmax)
    {
        CALL_AND_RETHROW(initialise_common());
        CALL_AND_RETHROW(initialise_indexing(wmax));
    }

    void initialise_cutoff(const T& wmax, const T& maxfactor)
    {
        CALL_AND_RETHROW(initialise_common());
        CALL_AND_RETHROW(initialise_indexing(wmax, maxfactor));
    }

    void initialise_cutoff(const T& wmax, const std::vector<T>& maxfactor, const std::vector<T>& scalefactor)
    {
        CALL_AND_RETHROW(initialise_common());
        CALL_AND_RETHROW(initialise_indexing(wmax, maxfactor, scalefactor));
    }

    template <typename terms_type>
    void initialise(const terms_type& o, const size_t& L)
    {
        CALL_AND_RETHROW(initialise_terms(o));
        CALL_AND_RETHROW(initialise_common());
        CALL_AND_RETHROW(initialise_indexing(L));
    }
    
    template <typename terms_type>
    void initialise(const terms_type& o, const T& wmax)
    {
        CALL_AND_RETHROW(initialise_terms(o));
        CALL_AND_RETHROW(initialise_common());
        CALL_AND_RETHROW(initialise_indexing(wmax));
    }

    template <typename terms_type>
    void initialise(const terms_type& o, const T& wmax, const T& maxfactor)
    {
        CALL_AND_RETHROW(initialise_terms(o));
        CALL_AND_RETHROW(initialise_common());
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
            std::cerr << nuk[i] << std::endl;
        }
    
        m_indexing.clear();
        CALL_AND_HANDLE(m_indexing.initialise(nuk, wmax), "Failed to initialise heom operator object.  Failed to initialise the basis indexing object.");
        m_indexing_initialised = true;
    }


    void initialise_indexing(const T& wmax, const std::vector<T>& maxfactor, const std::vector<T>& scalefactor)
    {
        std::vector<T> nuk(m_K);
        size_t counter = 0;
        T numin = -1;

        for(size_t i = 0; i < m_Kb.size(); ++i)
        {
            for(size_t k=0; k < m_Kb[i]; ++k)
            {
                T nui = linalg::real(m_nuk[counter]);
                if(numin > nui || numin < 0){numin=nui;}
            }
        }

        for(size_t i = 0; i < m_Kb.size(); ++i)
        {
            for(size_t k=0; k < m_Kb[i]; ++k)
            {
                nuk[counter] = linalg::real(m_nuk[counter])*scalefactor[i];
                if(nuk[counter] > maxfactor[i]*nuk[0]){nuk[counter] = maxfactor[i]*numin;}
                std::cerr << nuk[counter] << std::endl;
                ++counter;
            }
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

        m_jn.clear();
        m_terms.clear();
        m_Hs.clear();
        m_S.clear();
        m_Srho.clear();
        m_rhoS.clear();
    }



    void set_interactions(const linalg::matrix<complex_type> H, const linalg::matrix<complex_type>& S)
    {
        try
        {
            m_S.resize(1);  m_S[0] = S;
            m_Hs = H;
    
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to initialise sparse heom object.");
        }
    }

    void set_interactions(const linalg::matrix<complex_type> H, const std::vector<linalg::matrix<complex_type>>& S)
    {
        ASSERT(S.size() == m_terms.size(), "The number of bath correlation terms and system contributions to the system-bath coupling operators is not consistent.");
        size_t K = 0;
        m_Hs = H;

        //first determine the total number of modes
        for(size_t i=0; i < m_terms.size(); ++i)
        {
            for(size_t ti = 0; ti < m_terms[i].size(); ++ti)
            {
                K += m_terms[i][ti].nmodes();
                ASSERT(m_terms[i][ti].has_diagonal_gamma(), "sparse heom currently only support standard HEOM operators.");
            }
        }

        m_S = S;

        //now set up the frequency array and the two coupling operator terms
        m_jn.resize(K);
        m_c1k.resize(K);    m_c1k.fill_zeros();
        m_c2k.resize(K);    m_c2k.fill_zeros();
        m_c3k.resize(K);    m_c2k.fill_zeros();
        size_t counter = 0;
        m_K = K;
        m_nhilb = H.shape(0);

        for(size_t i=0; i < m_terms.size(); ++i)
        {
            for(size_t ti = 0; ti < m_terms[i].size(); ++ti)
            {
                for(size_t j=0; j < m_terms[i][ti].nmodes(); ++j)
                {
                    T prefactor = std::sqrt(linalg::abs(linalg::real(m_terms[i][ti].sk(j))));

                    if(m_terms[i][ti].has_symmetric_contribution(j))
                    {
                        m_c1k[counter] += (m_terms[i][ti].sk(j)/prefactor);
                        m_c2k[counter] += (m_terms[i][ti].sk(j)/prefactor);
                    }

                    if(m_terms[i][ti].has_asymmetric_contribution(j))
                    {                
                        m_c1k[counter] += ( complex_type(0, 1)*m_terms[i][ti].ak(j)/prefactor);
                        m_c2k[counter] -= ( complex_type(0, 1)*m_terms[i][ti].ak(j)/prefactor);
                    }

                    m_c3k[counter] += prefactor;
                    ++counter;
                }
            }
        }

        m_indexing_initialised = false;
        m_interactions_initialised = true;
    }


    void add_bath(const heom_bath<T>& bath)
    {   
        m_terms.push_back(bath);
        m_indexing_initialised = false;
    }

    size_t nados() const{return m_indexing.nstates();}
    size_t nmodes() const{return m_K;}

public:
    /*
     * Function for constructing the general HEOM evolution operator for a system interacting with a set of
     * baths with differing system coupling operators.  Here the user must specify as many system operators
     * as there are distinct terms stored by this object.
     */

    void apply(const linalg::vector<complex_type>& _A, linalg::vector<complex_type>& _B, complex_type coeff = complex_type(1, 0))
    {
        try
        {
            ASSERT(_A.size() == _B.size(), "Invalid input array sizes.");

            size_t nhilb = m_nhilb;
            size_t _nados = nados();
            _B.fill_zeros();

            auto A = _A.reinterpret_shape(_nados, nhilb, nhilb);
            auto B = _B.reinterpret_shape(_nados, nhilb, nhilb);


            m_jn.resize(m_K);
            m_Srho.resize(nhilb, nhilb);    m_Srho.fill_zeros();
            m_rhoS.resize(nhilb, nhilb);    m_rhoS.fill_zeros();

            for(size_t iado=0; iado < _nados; ++iado)
            {
                //we want to apply all operators that act on the Ai element in each step.  This minimises the number of times
                //that we need to evaluate the action of Srho and rhoS
                auto Ai = A[iado];

                m_indexing.state_at(iado, m_jn);
                complex_type gamma_n = complex_type(0, -1)*diagonal_bath_term(m_jn);

                //apply the diagonal term
                {
                    B[iado] += gamma_n * Ai;
                    m_rhoS = m_Hs*Ai;
                    B[iado] += m_rhoS;     
                    m_rhoS = Ai*m_Hs;
                    B[iado] -= m_rhoS;     
                }

                //now apply the raising and lowering terms
                size_t k = 0;
                for(size_t i=0; i < m_terms.size(); ++i)
                {
                    bool rhoS_evaluated = false;

                    for(size_t ti = 0; ti < m_terms[i].size(); ++ti)
                    {
                        for(size_t j=0; j < m_terms[i][ti].nmodes(); ++j)
                        {
                            if(m_indexing._contains_lowered_state(iado, k))
                            {
                                size_t ind = m_indexing._get_lowered_index(iado, k);
                                if(!rhoS_evaluated)
                                {
                                    m_Srho = m_S[i]*Ai;
                                    m_rhoS = Ai*m_S[i];
                                    rhoS_evaluated = true;
                                }
                
                                B[ind] += (sqrt(m_jn[k]*1.0))*m_c3k[k]*(m_Srho - m_rhoS);
                            }

                            if(m_indexing._contains_raised_state(iado, k))
                            {
                                size_t ind = m_indexing._get_raised_index(iado, k);
                                if(!rhoS_evaluated)
                                {
                                    m_Srho = m_S[i]*Ai;
                                    m_rhoS = Ai*m_S[i];
                                    rhoS_evaluated = true;
                                }
                
                                B[ind] += (sqrt(m_jn[k]+1.0))*(m_c1k[k]*m_Srho - m_c2k[k]*m_rhoS);
                            }
                            

                            ++k;
                        }
                        
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

    void initialise_terms(const heom_bath<T>& o)
    {
        try
        {
            m_terms.resize(1);
            m_terms[0] = o;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to initialise terms array.");
        }
    }


    void initialise_terms(const std::vector<heom_bath<T>>& o)
    {
        try
        {
            m_terms = o;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to initialise terms array.");
        }
    }

    void initialise_common()
    {
        m_K = 0;
        m_Kb.resize(m_terms.size());
        m_diagonal_bath = true;
        m_has_truncation_contribution = false;
        m_has_constant_truncation_contribution = false;
        m_has_non_constant_truncation_contribution = false;

        for(size_t i=0; i < m_terms.size(); ++i)
        {
            m_Kb[i] = 0;
            for(size_t ti = 0; ti < m_terms[i].size(); ++ti)
            {
                m_K += m_terms[i][ti].nmodes();
                m_Kb[i] += m_terms[i][ti].nmodes();
                m_diagonal_bath = m_diagonal_bath && m_terms[i][ti].has_diagonal_gamma();

                for(size_t j=0; j < m_terms[i][ti].nmodes(); ++j)
                {
                    m_has_asymmetric_coupling = m_has_asymmetric_coupling || m_terms[i][ti].has_asymmetric_contribution(j);
                }
            }
            if(m_terms[i].has_truncation())
            {
                m_has_truncation_contribution = true;
                if(m_terms[i].has_constant_truncation())
                {
                    m_has_constant_truncation_contribution = true;
                }
                else
                {
                    m_has_non_constant_truncation_contribution = true;
                }
            }
        }

        m_nuk.resize(m_K);
        size_t counter = 0;
        for(size_t i=0; i < m_terms.size(); ++i)
        {
            for(size_t ti = 0; ti < m_terms[i].size(); ++ti)
            {
                for(size_t j=0; j < m_terms[i][ti].nmodes(); ++j)
                {
                    m_nuk[counter] = m_terms[i][ti].nuk(j);
                    ++counter;
                }
            }
        }
    }

};

#endif

