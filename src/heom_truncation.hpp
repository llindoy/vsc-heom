#ifndef HEOM_TRUNCATION_HPP
#define HEOM_TRUNCATION_HPP

#include "bath_correlation_terms.hpp"
#include "superoperator_utilities.hpp"

template <typename T>
class heom_truncation
{
public:
    using complex_type = linalg::complex<T>;

    heom_truncation() {}
    virtual ~heom_truncation(){}
    virtual void initialise(std::shared_ptr<superoperator::superoperator<complex_type>> /* Ls */){}

    virtual void apply(std::shared_ptr<superoperator::superoperator<complex_type>> /* S */, const truncated_occupation_number_basis& /* index */, const std::vector<complex_type>& /* nuk */, const std::vector<size_t>& /* nt */, size_t /* iado */, linalg::matrix<complex_type>& /* M */) {}
    virtual bool constant_truncation() const {return true;}

    virtual void apply(const superoperator::diagonal_operator_superoperator<complex_type>& S, linalg::matrix<complex_type>& M ){}
};


template <typename T> 
class standard_truncation : public heom_truncation<T>
{
public:
    using complex_type = linalg::complex<T>;

    standard_truncation() : heom_truncation<T>() {}
    standard_truncation(const T& delta) : heom_truncation<T>(), m_delta(delta) {}
    ~standard_truncation(){}

    const T& delta() const{return m_delta;}
    T& delta() {return m_delta;}

    void apply(std::shared_ptr<superoperator::superoperator<complex_type>> S , const truncated_occupation_number_basis& /* index */, const std::vector<complex_type>& /* nuk */, const std::vector<size_t>& /* nt */, size_t /* iado */, linalg::matrix<complex_type>& M ) final
    {
        try
        {
            linalg::csr_matrix<complex_type> Sm;
            linalg::matrix<complex_type> So;
            S->commutator(Sm);
            superoperator::csr_to_dense(Sm, So);

            M += complex_type(0, -1)*m_delta*So*So;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct system operator matrix.");
        }
    }

    void apply(const superoperator::diagonal_operator_superoperator<complex_type>& S, linalg::matrix<complex_type>& M ) final
    {
        try
        {
            linalg::csr_matrix<complex_type> Sm;
            linalg::matrix<complex_type> So;
            S.commutator(Sm);
            superoperator::csr_to_dense(Sm, So);

            M += complex_type(0, -1)*m_delta*So*So;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct system operator matrix.");
        }
    }
protected:
    T m_delta;
};

template <typename T> 
class nakajima_zwanzig_truncation : public heom_truncation<T>
{
public:
    using complex_type = linalg::complex<T>;

    nakajima_zwanzig_truncation() : heom_truncation<T>(){}
    nakajima_zwanzig_truncation(const correlation_terms<T>& pade, const correlation_terms<T>& mats) : heom_truncation<T>(), m_pade(pade), m_mats(mats){}
    nakajima_zwanzig_truncation(correlation_terms<T>&& pade, correlation_terms<T>&& mats) : heom_truncation<T>(), m_pade(std::forward(pade)), m_mats(std::forward(mats)){}

    ~nakajima_zwanzig_truncation(){}

    void set_pade(const correlation_terms<T>& pade){m_pade = pade;}
    void set_pade(correlation_terms<T>&& pade){m_pade = std::move(pade);}
    void set_mats(const correlation_terms<T>& mats){m_mats = mats;}
    void set_mats(correlation_terms<T>&& mats){m_mats = std::move(mats);}

    const bool& apply_terminator() const{return m_apply_terminator;}
    bool& apply_terminator(){return m_apply_terminator;}

    virtual void initialise(std::shared_ptr<superoperator::superoperator<complex_type>> H )
    {
        ASSERT(m_pade.has_diagonal_gamma() && m_mats.has_diagonal_gamma(), "Currently doesn't support non-diagonal gamma.");
        linalg::csr_matrix<complex_type> L;  
        H->commutator(L);
        linalg::matrix<complex_type> Ld;    superoperator::csr_to_dense(L, Ld);
        linalg::matrix<complex_type> mLs;   

        linalg::eigensolver<linalg::hermitian_matrix<complex_type>> solver(Ld);
        solver(Ld, m_E, m_U);

        m_initialised = true;
        m_diag.resize(m_E.size(0), m_E.size(1));
        m_temp.resize(m_E.size(0), m_E.size(1));
        m_temp2.resize(m_E.size(0), m_E.size(1));

        m_Lkm.resize(m_E.size(0), m_E.size(1));
        m_Lkp.resize(m_E.size(0), m_E.size(1));
        m_Cum.resize(m_E.size(0), m_E.size(1));
    }

    void apply(std::shared_ptr<superoperator::superoperator<complex_type>> S , const truncated_occupation_number_basis& index, const std::vector<complex_type>& nuk, const std::vector<size_t>& n, size_t iado, linalg::matrix<complex_type>& M) final
    {
        try
        {
            ASSERT(m_initialised, "Have not initialised the NZ truncation object.");

            //evaluate the gamma term for the current ado
            complex_type gamma(0, 0);
            for(size_t k = 0; k < nuk.size(); ++k){gamma += nuk[k] * static_cast<T>(n[k]);}


            //first we construct Lkm as the commutator
            S->commutator(m_Sm);
            S->anticommutator(m_Sp);
            superoperator::csr_to_dense(m_Sm, m_Lkm);
        
            m_Cum.fill_zeros();
            //deal with the truncation of the matsubara series expansion
            //first we add on the contributions from a high order pade approximant in order to accurately approximate the full nz series
            for(size_t k=0; k < m_pade.nmodes(); ++k)
            {
                //first we construct the Lkp term as that depends on k
                form_lkp(m_pade, k);

                //now we can form the matrix inverse term

                for(size_t i = 0; i < m_E.size(0); ++i)
                {
                    m_diag(i, i) = 1.0/(gamma + m_pade.nuk(k) + complex_type(0, m_E(i, i)));
                }

                m_temp = m_U*m_diag;
                m_temp2 = m_temp*linalg::adjoint(m_U);
                m_temp = m_temp2*m_Lkp;
                m_Cum += m_temp;


            }

            //now subtract off the included matsubara terms
            for(size_t k=0; k < m_mats.nmodes(); ++k)
            {
                //first we construct the Lkp term as that depends on k
                form_lkp(m_mats, k);

                //now we can form the matrix inverse term

                for(size_t i = 0; i < m_E.size(0); ++i)
                {
                    m_diag(i, i) = 1.0/(gamma + m_mats.nuk(k) + complex_type(0, m_E(i, i)));
                }

                m_temp = m_U*m_diag;
                m_temp2 = m_temp*linalg::adjoint(m_U);
                m_temp = m_temp2*m_Lkp;
                m_Cum -= m_temp;

                //now we attempt to apply the terminator
                if(m_apply_terminator)
                {

                }
            }

            M += m_Lkm*m_Cum; 
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct system operator matrix.");
        }
    }
    bool constant_truncation() const final {return false;}

protected:
    void form_lkp(const correlation_terms<T>& corr, size_t k)
    {
        auto rp = m_Sm.rowptr();
        auto ci = m_Sm.colind();
        auto smb = m_Sm.buffer();
        auto spb = m_Sp.buffer();

        for(size_t i = 0; i < m_Lkp.size(0); ++i)
        {
            for(size_t c = rp[i]; c < static_cast<size_t>(rp[i+1]); ++c)
            {
                m_Lkp(i, ci[c]) = ((corr.has_symmetric_contribution(k) ? corr.sk(k) : 0.0)*smb[c] + (corr.has_asymmetric_contribution(k) ? complex_type(0, 1)*corr.ak(k) : 0.0)*spb[c]);
            }
        }
    }

protected:
    linalg::matrix<complex_type> m_U;
    linalg::matrix<complex_type> m_temp;
    linalg::matrix<complex_type> m_temp2;
    linalg::diagonal_matrix<T> m_E;
    linalg::diagonal_matrix<complex_type> m_diag;
    linalg::csr_matrix<complex_type> m_Sm;
    linalg::csr_matrix<complex_type> m_Sp;
    linalg::matrix<complex_type> m_Lkp;
    linalg::matrix<complex_type> m_Lkm;
    linalg::matrix<complex_type> m_Cum;

    correlation_terms<T> m_pade;
    correlation_terms<T> m_mats;
    bool m_initialised;
    bool m_apply_terminator;
};


#endif  //HEOM_TRUNCATION_HPP//

