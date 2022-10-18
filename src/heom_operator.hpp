#ifndef HEOM_OPERATOR_HPP
#define HEOM_OPERATOR_HPP

#include <linalg/linalg.hpp>

#include "occupation_number_basis_indexing.hpp"
#include "heom_bath_terms.hpp"
#include "superoperator_utilities.hpp"


//an object for constructing generic HEOM bath operators.  For an arbitrary spectral density the generalised heom can be written in a system bath form 
//with an operator acting on the bath degrees of freedom (no action on the system degrees of freedom) and up to two terms that couple the system and bath degrees
//of freedom.  A term accounting for the symmetric part of the bath correlation function and potentially a term accounting for the asymmetric part.  
template <typename T>
class heom_bath_operator
{
public:
    using complex_type = linalg::complex<T>;

protected:
    truncated_occupation_number_basis m_indexing;

    size_t m_K;
    std::vector<size_t> m_Kb;
    bool m_diagonal_bath;
    bool m_has_asymmetric_coupling;

    bool m_has_truncation_contribution;
    bool m_has_non_constant_truncation_contribution;
    bool m_has_constant_truncation_contribution;

    bool m_indexing_initialised;
    bool m_rescale_couplings;
    std::vector<complex_type> m_nuk;

    std::vector<heom_bath<T>> m_terms;
public: 
    heom_bath_operator() : m_K(0), m_diagonal_bath(false), m_has_asymmetric_coupling(false), m_indexing_initialised(false), m_rescale_couplings(true) {}

    heom_bath_operator(const heom_bath<T>& o, const size_t& L) : heom_bath_operator()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }
    heom_bath_operator(const std::vector<heom_bath<T>>& o, const size_t& L) : heom_bath_operator()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }
    heom_bath_operator(const std::vector<std::vector<heom_bath<T>>>& o, const size_t& L) : heom_bath_operator()
    {
        CALL_AND_HANDLE(initialise(o, L), "Failed to construct heom operator object.  Failed to initialise bath correlation function.");
    }

    heom_bath_operator(const heom_bath_operator& o) = default;
    heom_bath_operator(heom_bath_operator&& o) = default;

    heom_bath_operator& operator=(const heom_bath_operator& o) = default;
    heom_bath_operator& operator=(heom_bath_operator&& o) = default;

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
        m_terms.clear();
    }


    void add_bath(const heom_bath<T>& bath)
    {   
        m_terms.push_back(bath);
        m_indexing_initialised = false;
    }

    size_t nados() const{return m_indexing.nstates();}
    size_t nmodes() const{return m_K;}
    size_t nbaths() const{return m_terms.size();}

public:
    /*
     *  Functions for creating the HEOM system operator.  This is simply the Liouvillian plus the markovian correction terms
     */
    void system_operator(std::shared_ptr<superoperator::superoperator<complex_type>> H, const std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>>& S, linalg::csr_matrix<complex_type>& L, T prune_tol = -1) const
    {
        try
        {
            if(m_has_constant_truncation_contribution)
            {
                linalg::matrix<complex_type> M;
                CALL_AND_RETHROW(system_operator(H, S, M));
                superoperator::dense_to_csr(M, L, prune_tol);
            }
            else
            {
                H->commutator(L);
                if(prune_tol > 0)
                {
                    L.prune(prune_tol);
                }
            }

            //then we make the csr matrix for the system variables be dense
            if(m_has_non_constant_truncation_contribution)
            {
                linalg::matrix<complex_type> M;
                superoperator::csr_to_dense(L, M);
                superoperator::dense_to_csr(M, L, -1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct system operator matrix.");
        }
    }

    void system_operator(std::shared_ptr<superoperator::superoperator<complex_type>> H, std::shared_ptr<superoperator::superoperator<complex_type>> S, linalg::csr_matrix<complex_type>& L, T prune_tol = -1) const
    {
        std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>> Sv(1);
        Sv[0] = S;
        CALL_AND_RETHROW(system_operator(H, Sv, L, prune_tol));
    }

    void system_operator(std::shared_ptr<superoperator::superoperator<complex_type>> H, const std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>>& S, linalg::matrix<complex_type>& M) const
    {
        try
        {
            ASSERT(S.size() == m_terms.size(), "The number of coupling operators must be the same as the number of bath objects.");
            

            linalg::csr_matrix<complex_type> L;
            H->commutator(L);

            linalg::matrix<complex_type> So;
            superoperator::csr_to_dense(L, M);

            std::vector<size_t> nt(m_K, 0);
            for(size_t tt = 0; tt < S.size(); ++tt)
            {
                if(m_terms[tt].has_constant_truncation())
                {
                    m_terms[tt].truncation()->apply(S[tt], m_indexing, m_nuk, nt, 0, M);
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct system operator matrix.");
        }
    }

    void system_operator(std::shared_ptr<superoperator::superoperator<complex_type>> H, std::shared_ptr<superoperator::superoperator<complex_type>> S, linalg::matrix<complex_type>& L) const
    {
        std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>> Sv(1);
        Sv[0] = S;
        CALL_AND_RETHROW(system_operator(H, Sv, L));
    }
    
public:
    /*
     * Function for constructing the general HEOM evolution operator for a system interacting with a set of
     * baths with differing system coupling operators.  Here the user must specify as many system operators
     * as there are distinct terms stored by this object.
     */
    void operator()(std::shared_ptr<superoperator::superoperator<complex_type>> H, const std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>>& S, linalg::csr_matrix<complex_type>& op, T prune_tol = -1)
    {
        try
        {
            ASSERT(S.size() == m_terms.size(), "The number of coupling operators must be the same as the number of bath objects.");

            //if we have a truncation contribution then we initialise the truncation objects using the system liouvillian
            if(m_has_truncation_contribution)
            {
                for(size_t i = 0; i < m_terms.size(); ++i)
                {
                    if(m_terms[i].has_truncation())
                    {
                        m_terms[i].truncation()->initialise(H);
                    }
                }
            }

            //build all superoperators required for these calculations
            linalg::csr_matrix<complex_type> L;  
            system_operator(H, S, L, prune_tol);

            size_t nhilb = H->hilb_dim();
            size_t nliouv = nhilb*nhilb;

            //store whether or not the liouvillian contains the diagonal term on this row.
            std::vector<bool> L_has_diagonal(nliouv);
            for(size_t iliouv = 0; iliouv < nliouv; ++iliouv)
            {
                L_has_diagonal[iliouv] =  L.contains_diagonal(iliouv);
            }
    
            std::vector<linalg::csr_matrix<complex_type>> Sm(S.size());
            std::vector<linalg::csr_matrix<complex_type>> Sp(S.size());

            //here we don't exploit the fact that some of the commutator terms are zero in the sparsity structure of the evolution operator, as this complicates the indexing significantly
            for(size_t i = 0; i < S.size(); ++i)
            {
                ASSERT(S[i]->hilb_dim() == nhilb, "The coupling operators do not have the same dimension as the Hamiltonian.");

                S[i]->commutator(Sm[i]);
                S[i]->anticommutator(Sp[i]);

                ASSERT(Sm[i].topology() == Sp[i].topology(), "The topology of the commutator and anticommutator superoperators are not the same.");
            }

            size_t _nados = nados();
            std::vector<size_t> nt(m_K, 0);

            std::vector<std::vector<size_t>> nkskip(m_terms.size());
            size_t skip_counter = 0;
            for(size_t i = 0; i < m_terms.size(); ++i)
            {
                nkskip[i].resize(m_terms[i].size(), skip_counter);
                for(size_t j = 0; j < m_terms[i].size(); ++j)
                {
                    nkskip[i][j] = skip_counter;
                    skip_counter += m_terms[i][j].nmodes();
                }
            }

            //now that we have formed each of the superoperators we can actually go through and construct the HEOM evolution matrix
            //as a csr matrix.  Here we will make the ado index be the slow index and the Liouville space index the fast index.
            //For each term we add in the terms that couple to lowered states followed by the lower-triangle bath-bath coupling terms.  
            //This is followed by the lower triangular Liouvillian terms, then the diagonal terms (both system and bath), the upper 
            //triangular liouvillian terms, the upper-triangular bath-bath terms and finally the terms coupling to raised states.
            CALL_AND_HANDLE(op.resize(_nados*nliouv, _nados*nliouv), "Failed to resize HEOM evolution operator matrix.");
            
            auto rowptr = op.rowptr();  rowptr[0] = 0;

            //first we work out the total number of non-zero elements in the evolution operator and set up the rowptr array.
            size_t nnz = 0;
            
            for(size_t iado=0; iado < _nados; ++iado)
            {
                for(size_t iliouv = 0; iliouv < nliouv; ++iliouv)
                {
                    //add in the number of non-zeros that account for the raising and lowering terms.   These terms are either coupled to
                    //the commutator or anticommutator which here have the same sparsity pattern.   
                    for(size_t tt = 0; tt < m_terms.size(); ++tt)
                    {
                        for(size_t ti = 0; ti < m_terms[tt].size(); ++ti)
                        {
                            size_t Ki = m_terms[tt][ti].nmodes();
                            for(size_t k1 = 0; k1 < Ki; ++k1)
                            {
                                size_t _k1 = k1+nkskip[tt][ti];
                                //add in the term connected to increased states that is always present
                                nnz += (m_indexing.contains_raised_state(iado, _k1) ? Sm[tt].nnz_in_row(iliouv) : 0);

                                //attempt to add in the term connected to the lowered state.  This is added if we have either symmetric or asymmetic terms
                                if(m_terms[tt][ti].has_symmetric_contribution(k1) || m_terms[tt][ti].has_symmetric_contribution(k1))
                                {
                                    nnz +=  (m_indexing.contains_lowered_state(iado, _k1) ? Sm[tt].nnz_in_row(iliouv) : 0);
                                }
                            }
                        }   
                    }

                    //now add in all of the bath states that arise from adjacent terms.  These terms are diagonal in the system operator so we only 
                    //add a single term per row per contribution.
                    for(size_t tt = 0; tt < m_terms.size(); ++tt)
                    {
                        for(size_t ti = 0; ti < m_terms[tt].size(); ++ti)
                        {
                            if(!m_terms[tt][ti].has_diagonal_gamma())
                            {
                                size_t Ki = m_terms[tt][ti].nmodes();
                                for(size_t k1 = 0; k1 < Ki; ++k1)
                                {
                                    for(size_t k2 = 0; k2 < Ki; ++k2)
                                    {
                                        if(k1 != k2 && m_terms[tt][ti].nonzero_gamma(k1, k2))
                                        {
                                            nnz += (m_indexing.contains_adjacent_state(iado, k1+nkskip[tt][ti], k2+nkskip[tt][ti])  ? 1 : 0);
                                        }
                                    }
                                }
                            }
                        }   
                    }

                    //now we deal with the diagonal (in ado space contributions).  If the liouvillian does not contain the diagonal element then 
                    //we need to add in the purely diagonal element
                    nnz += L.nnz_in_row(iliouv) + (L_has_diagonal[iliouv] ? 0 : 1);

                    //first add on the lowered index terms
                    rowptr[iado*nliouv+iliouv+1] = nnz;
                }
            }
            std::cerr << nnz << std::endl;
            CALL_AND_HANDLE(op.resize(nnz), "Failed to resize HEOM evolution operator matrix.");
        
            //now we actually go ahead and build the operator.
            auto colind = op.colind();
            auto buffer = op.buffer();

            //now we actually go ahead and set up the correct matrix buffers.
            size_t counter = 0; 
            for(size_t iado=0; iado < _nados; ++iado)
            {
                m_indexing.state_at(iado, nt);

                //now set up the diagonal term for this ado.  This includes any non-constant truncation terms
                linalg::csr_matrix<complex_type> Ll(L);  
                if(m_has_non_constant_truncation_contribution)
                {
                    std::cerr << iado << " " << _nados << std::endl;
                    linalg::matrix<complex_type> M;
                    superoperator::csr_to_dense(Ll, M);

                    for(size_t tt = 0; tt < S.size(); ++tt)
                    {
                        if(m_terms[tt].has_truncation() && !m_terms[tt].has_constant_truncation())
                        {
                            m_terms[tt].truncation()->apply(S[tt], m_indexing, m_nuk, nt, iado, M);
                        }
                    }
                    superoperator::dense_to_csr(M, Ll, -1.0);
                }
            
                complex_type gamma_n = diagonal_bath_term(nt);

                //add in the terms that couple the current ado to ados with lower the number of excitations
                for(size_t iliouv = 0; iliouv < nliouv; ++iliouv)
                {
                    //first add the terms where the column index is smaller than i.  The lowering terms.
                    for(size_t tt = 0; tt < m_terms.size(); ++tt)
                    {
                        for(size_t ti = 0; ti < m_terms[tt].size(); ++ti)
                        {
                            size_t Ki = m_terms[tt][ti].nmodes();
                            for(size_t k1 = 0; k1 < Ki; ++k1)
                            {
                                size_t _k1 = k1+nkskip[tt][ti];
                                if(m_indexing.contains_lowered_state(iado, _k1))
                                {
                                    size_t ind = m_indexing.get_lowered_index(iado, _k1);
    
                                    T prefactor = 1.0;
                                    if(m_rescale_couplings && m_terms[tt][ti].has_symmetric_contribution(k1))   
                                    {
                                        prefactor = std::sqrt(linalg::abs(linalg::real(m_terms[tt][ti].sk(k1))));
                                    }
                                
                                    for(size_t rpi = static_cast<size_t>(Sm[tt].rowptr()[iliouv]); rpi < static_cast<size_t>(Sm[tt].rowptr()[iliouv+1]); ++rpi)
                                    {   
                                        complex_type val(0, 0);
                                        bool term_added =false;
                                        if(m_terms[tt][ti].has_symmetric_contribution(k1))
                                        {
                                            //val += m_terms[tt][ti].sk(k1)/prefactor*std::pow(static_cast<T>(nt[_k1]), 0.05)*Sm[tt].buffer()[rpi];
                                            val += m_terms[tt][ti].sk(k1)/prefactor*std::sqrt(static_cast<T>(nt[_k1]))*Sm[tt].buffer()[rpi];
                                            term_added = true;
                                        }
                                        if(m_terms[tt][ti].has_asymmetric_contribution(k1))
                                        {
                                            //val += complex_type(0, 1)*m_terms[tt][ti].ak(k1)/prefactor*std::pow(static_cast<T>(nt[_k1]), 0.05)*Sp[tt].buffer()[rpi];
                                            val += complex_type(0, 1)*m_terms[tt][ti].ak(k1)/prefactor*std::sqrt(static_cast<T>(nt[_k1]))*Sp[tt].buffer()[rpi];
                                            term_added = true;
                                        }
                
                                        if(term_added)
                                        {
                                            buffer[counter] = val;
                                            colind[counter] = ind*nliouv + Sm[tt].colind()[rpi];
                                            ++counter;
                                        }
                                    }
                                }
                            }
                        }   
                    }

                    //add the terms that preserve the total number of excitations but where the column index is smaller than the row index
                    for(size_t tt = 0; tt < m_terms.size(); ++tt)
                    {
                        for(size_t ti = 0; ti < m_terms[tt].size(); ++ti)
                        {
                            if(!m_terms[tt][ti].has_diagonal_gamma())
                            {
                                size_t Ki = m_terms[tt][ti].nmodes();
                                for(size_t k1 = 0; k1 < Ki; ++k1)
                                {
                                    size_t _k1 = k1+nkskip[tt][ti];
                                    for(size_t k2 = Ki-1; k2 >= k1+1; --k2)
                                    {
                                        size_t _k2 = k2+nkskip[tt][ti];
                                        if(m_terms[tt][ti].nonzero_gamma(k1, k2))
                                        {
                                            if(m_indexing.contains_adjacent_state(iado, _k1, _k2))
                                            {
                                                size_t ind = m_indexing.get_adjacent_index(iado, _k1, _k2);

                                                buffer[counter] = complex_type(0, -1)*m_terms[tt][ti].gamma(k1, k2)*std::sqrt(nt[_k1]*(nt[_k2]+1.0));
                                                colind[counter] = ind*nliouv+iliouv;
                                                ++counter;
                                            }
                                        }
                                    }
                                }
                            }
                        }   
                    }

                    {
                            


                        //now we add in the terms that are diagonal in ado space
                        size_t rpi = 0;
                        for(rpi = static_cast<size_t>(Ll.rowptr()[iliouv]); rpi < static_cast<size_t>(Ll.rowptr()[iliouv+1]) && static_cast<size_t>(Ll.colind()[rpi]) < iliouv; ++rpi)
                        {
                            buffer[counter] = Ll.buffer()[rpi];
                            colind[counter] = iado*nliouv + Ll.colind()[rpi];
                            ++counter;
                        }
                        
                        //then add the diagonal terms
                        {
                            
                            buffer[counter] = complex_type(0, -1)*gamma_n;
                            if(static_cast<size_t>(Ll.colind()[rpi]) == iliouv)
                            {
                                buffer[counter] += Ll.buffer()[rpi];
                                ++rpi;
                            }
                            colind[counter] = iado*nliouv+iliouv;
                            ++counter;
                        }
                        for(; rpi < static_cast<size_t>(Ll.rowptr()[iliouv+1]); ++rpi)
                        {
                            buffer[counter] = Ll.buffer()[rpi];
                            colind[counter] = iado*nliouv + Ll.colind()[rpi];
                            ++counter;
                        }
                    }

                    //now we the excitation preserving terms with column index larger than the row index
                    for(size_t tt = 0; tt < m_terms.size(); ++tt)
                    {
                        for(size_t ti = 0; ti < m_terms[tt].size(); ++ti)
                        {
                            if(!m_terms[tt][ti].has_diagonal_gamma())
                            {
                                size_t Ki = m_terms[tt][ti].nmodes();
                                for(size_t j2 = 0; j2 < Ki; ++j2)
                                {
                                    size_t k2 = Ki - (j2+1);
                                    size_t _k2 = k2+nkskip[tt][ti];
                                    for(size_t k1 = k2+1; k1 < Ki; ++k1)
                                    {
                                        size_t _k1 = k1+nkskip[tt][ti];
                                        if(m_terms[tt][ti].nonzero_gamma(k1, k2))
                                        {
                                            if(m_indexing.contains_adjacent_state(iado, _k1, _k2))
                                            {
                                                size_t ind = m_indexing.get_adjacent_index(iado, _k1, _k2);

                                                buffer[counter] = complex_type(0, -1)*m_terms[tt][ti].gamma(k1, k2)*std::sqrt(nt[_k1]*(nt[_k2]+1.0));
                                                colind[counter] = ind*nliouv+iliouv;
                                                ++counter;
                                            }
                                        }
                                    }
                                }
                            }
                        }   
                    }

                    //now we add the terms that couple this ado to those with larger numbers of excitations
                    for(size_t tt = 0; tt < m_terms.size(); ++tt)
                    {
                        for(size_t ti = 0; ti < m_terms[tt].size(); ++ti)
                        {
                            size_t Ki = m_terms[tt][ti].nmodes();
                            for(size_t j1 = 0; j1 < Ki; ++j1)
                            {
                                size_t k1 = Ki - (j1+1);
                                size_t _k1 = k1+nkskip[tt][ti];
                                if(m_indexing.contains_raised_state(iado, _k1))
                                {
                                    size_t ind = m_indexing.get_raised_index(iado, _k1);
                                    
                                    T prefactor = 1;
                                    if(m_terms[tt][ti].has_symmetric_contribution(k1)  && m_rescale_couplings)
                                    {
                                        prefactor = std::sqrt(linalg::abs(linalg::real(m_terms[tt][ti].sk(k1))));
                                    }

                                    for(size_t rpi = static_cast<size_t>(Sm[tt].rowptr()[iliouv]); rpi < static_cast<size_t>(Sm[tt].rowptr()[iliouv+1]); ++rpi)
                                    {   
                                        if(m_terms[tt][ti].has_symmetric_contribution(k1))
                                        {
                                            //buffer[counter] = std::pow((nt[_k1]+1.0), 0.95)*prefactor*Sm[tt].buffer()[rpi];
                                            buffer[counter] = std::sqrt((nt[_k1]+1.0))*prefactor*Sm[tt].buffer()[rpi];
                                            colind[counter] = ind*nliouv + Sm[tt].colind()[rpi];
                                            ++counter;
                                        }
                                    }
                                }
                            }
                        }   
                    }
                }
            }
            if(prune_tol > 0)
            {
                op.prune(prune_tol);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct HEOM evolution operator.");
        }
    } 
    

    void operator()(std::shared_ptr<superoperator::superoperator<complex_type>> H, std::shared_ptr<superoperator::superoperator<complex_type>> S, linalg::csr_matrix<complex_type>& op, T prune_tol = -1)
    {
        std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>> Sv(1);
        Sv[0] = S;
        CALL_AND_RETHROW(this->operator()(H, Sv, op, prune_tol));
    }

protected:
    complex_type diagonal_bath_term(const std::vector<size_t>& n) const
    {
        complex_type ci(0, 0);
        for(size_t k = 0; k < m_K; ++k){ci += m_nuk[k] * static_cast<T>(n[k]);}
        return ci;
    }

protected:

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

