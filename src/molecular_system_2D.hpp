#ifndef MOLECULAR_BASIS_2D_HPP
#define MOLECULAR_BASIS_2D_HPP

#include "potential.hpp"
#include "molecular_system.hpp"
#include "hermite_dvr.hpp"


template <typename T> 
class molecular_system_2D
{
public:
    molecular_system_2D() : m_diagonalise_position(false){}
    molecular_system_2D(const rapidjson::Value& obj) : m_diagonalise_position(false){CALL_AND_RETHROW(load(obj));}

    template <typename ... Args> 
    molecular_system_2D(size_t nT, size_t nQ, Args&& ...args) : m_molsystem(std::forward<Args>(args)...), _nT(nT), _nQ(nQ) {}

    const size_t& Nm() const{return _nT;}
    size_t& Nm(){return _nT;}

    const size_t& NQ() const{return _nQ;}
    size_t& NQ(){return _nQ;}

    void construct_basis()
    {
        m_molsystem.construct_basis();
        //first set up the renormalised system Hamiltonian
        form_hamiltonian();
        linalg::matrix<T> Ut(m_Hp.shape(0), m_Hp.shape(1));  
        linalg::diagonal_matrix<T> Et(m_Hp.shape(0));  

        //and diagonalise it to find the relevant fbr states
        CALL_AND_HANDLE(_solver(m_Hp, Et, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

        _nT = 0;
        for(size_t i = 0; i < m_Hp.shape(0);++i)
        {
            T cmn1 = 1.0/219474.63;
            if(Et(i, i)/cmn1 < _Ep){++_nT;}
        }
        if(_nT < _nTmin){_nT = _nTmin;}
        m_E.resize(_nT);

        //now truncate the system to only include the nm states we wish to retain
        m_E.resize(_nT);  m_U.resize(_nT, m_Hp.size(0));
        for(size_t i = 0; i < _nT; ++i)
        {
            m_E(i) = Et(i, i);
            std::cerr << m_E(i)*27.2114 << std::endl;
            for(size_t j = 0; j < m_Hp.size(0); ++j)
            {
                m_U(i, j) = Ut(j, i);
            }
        }

        //now form the position operator in the polaritonic basis
        if(m_diagonalise_position)
        {
            //form the position operator in the polaritonic basis
            m_R.resize(_nT, _nT);   m_R.fill_zeros();
            size_t nM = m_molsystem.Nm();
            linalg::vector<T> xm = m_molsystem.x();
            diag.resize(nM*_nQ, nM*_nQ);
            for(size_t i = 0; i<nM; ++i)
            {
                for(size_t j = 0; j < _nQ; ++j)
                {
                    diag(i*_nQ+j, i*_nQ+j) = xm(i);                               
                }
            }
            transform_operator(diag, m_R);

            //and diagonalise it to find the relevant fbr states
            CALL_AND_HANDLE(_solver(m_R, m_Rv, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

            //now we need to tranform the transformation matrix and the polaritonic Hamiltonian
            {
                linalg::matrix<T> temp2 = linalg::trans(Ut)*m_U;
                m_U = temp2;
            }

            //now set up the Hamiltonian in the new basis
            {
                linalg::diagonal_matrix<T> Ev(m_E.size(),  m_E.size());
                for(size_t i =0; i < m_E.size(); ++i){Ev(i, i) = m_E(i);}
                linalg::matrix<T> temp2 = Ev*Ut;
                m_Hpf = linalg::trans(Ut)*temp2;
            }
        }
    }

    void hamiltonian_operator(linalg::matrix<T>& op, T /* renorm */ = 0)
    {
        op.resize(_nT, _nT);
        if(!m_diagonalise_position)
        {
            op.fill_zeros();
            for(size_t i =0; i < _nT; ++i){op(i, i) = m_E(i);}
        }
        else
        {
            for(size_t i = 0; i < _nT; ++i)
            {
                for(size_t j=0; j < _nT; ++j)
                {
                    op(i, j) = m_Hpf(i, j);
                }
            }
        }
    }

    const linalg::vector<T>& E() const{return m_E;}

    void R_operator(linalg::matrix<T>& op, T Rdisp=0)
    {
        op.resize(_nT, _nT);
        op.fill_zeros();

        if(!m_diagonalise_position)
        {
            size_t nM = m_molsystem.Nm();
            linalg::vector<T> xm = m_molsystem.x();
            diag.resize(nM*_nQ, nM*_nQ);
            for(size_t i = 0; i<nM; ++i)
            {
                for(size_t j = 0; j < _nQ; ++j)
                {
                    diag(i*_nQ+j, i*_nQ+j) = xm(i)-Rdisp;                               
                }
            }
            transform_operator(diag, op);
        }
        else
        {
            for(size_t i =0; i < _nT; ++i){op(i, i) = m_Rv(i)-Rdisp;}
        }
    }

    void R2_operator(linalg::matrix<T>& op)
    {
        op.resize(_nT, _nT);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::vector<T> xm = m_molsystem.x();
        diag.resize(nM*_nQ, nM*_nQ);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t j = 0; j < _nQ; ++j)
            {
                diag(i*_nQ+j, i*_nQ+j) = xm(i)*xm(i);                               
            }
        }
        transform_operator(diag, op);
    }

    void Q_operator(linalg::matrix<T>& op, T /* Rdisp */=0)
    {
        op.resize(_nT, _nT);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::vector<T> xc = m_herm.x();
        diag.resize(nM*_nQ, nM*_nQ);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t j = 0; j < _nQ; ++j)
            {
                diag(i*_nQ+j, i*_nQ+j) = xc(j)/std::sqrt(_w);
            }
        }
        transform_operator(diag, op);
    }

    void Q2_operator(linalg::matrix<T>& op)
    {
        op.resize(_nT, _nT);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::vector<T> xc = m_herm.x();
        diag.resize(nM*_nQ, nM*_nQ);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t j = 0; j < _nQ; ++j)
            {
                diag(i*_nQ+j, i*_nQ+j) = xc(j)*xc(j)/_w;                               
            }
        }
        transform_operator(diag, op);
    }

    void side_operator(linalg::matrix<T>& op)
    {
        op.resize(_nT, _nT);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::vector<T> xm = m_molsystem.x();
        diag.resize(nM*_nQ, nM*_nQ);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t j = 0; j < _nQ; ++j)
            {
                diag(i*_nQ+j, i*_nQ+j) = xm(i) < 0.0 ? 1.0 : 0.0;                               
            }
        }
        transform_operator(diag, op);

    }

    void rho0(linalg::matrix<T>& op, T beta, T renorm)
    {
        op.resize(_nT, _nT);
        op.fill_zeros();
        linalg::matrix<T> H, R2;
        hamiltonian_operator(H);
        R2_operator(R2);
        H = H+renorm*R2;

        temp.resize(_nT, _nT);
        linalg::matrix<T> Ut(_nT, _nT);  
        linalg::diagonal_matrix<T> Et(_nT);  

        //and diagonalise it to find the relevant fbr states
        CALL_AND_HANDLE(_solver(H, Et, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

        Et = elemental_exp(-beta*Et);

        temp = Ut*Et;
        op = temp*linalg::trans(Ut);
    }
protected:
    //here we form the cavity Hamiltonian in the hermite DVR basis and construct the mixed
    //molecule-system Hamiltonian
    void form_hamiltonian()
    {
        size_t nM = m_molsystem.Nm();
        
        m_herm.w0() = _w;


        m_Hp.resize(nM*_nQ, nM*_nQ);    m_Hp.fill_zeros();
        linalg::matrix<T> Hm;
        //first form the molecular Hamiltonian
        m_molsystem.hamiltonian_operator(Hm);
        linalg::vector<T> xm = m_molsystem.x();

        //now form the cavity harmonic oscillator Hamiltonian
        m_herm.Ndvr() = _nQ;
        m_herm.construct_quadrature();
        linalg::matrix<T> Hc(_nQ, _nQ);
        m_herm.H0(Hc);
        linalg::vector<T> xc = m_herm.x();

        //here we are working with frequency scaled positions for the cavity mode.  E.g. we are working with
        //x = sqrt(w_c) q.  So make sure we scale all of the coupling terms and self energy terms appropriately

        for(size_t i = 0; i < nM; ++i)
        {
            //first add on the purely molecular terms
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                for(size_t j = 0; j < _nQ; ++j)
                {
                    //add on the diagonal terms
                    size_t k = i*_nQ+j;
                    size_t k2 = i2*_nQ+j;
                    m_Hp(k, k2) += Hm(i,  i2);        //the matter terms

                }
            }

            //add on the dipole self energy terms
            for(size_t j = 0; j < _nQ; ++j)
            {   
                size_t k = i*_nQ+j;
                m_Hp(k, k) += 0.5*_c*_c/(_w*_w)*xm(i)*xm(i);        //the dipole self energy terms
            }


            //now we add on the purely cavity terms
            for(size_t j = 0; j < _nQ; ++j)
            {
                for(size_t j2 = 0; j2 < _nQ; ++j2)
                {
                    //add on the diagonal terms
                    size_t k = i*_nQ+j;
                    size_t k2 = i*_nQ+j2;
                    m_Hp(k, k2) += Hc(j,  j2);        //the matter terms - including the dipole self energy terms
                }
            }

            //now add in the molecule cavity coupling
            for(size_t j = 0; j < _nQ; ++j)
            {
                //add on the diagonal terms
                size_t k = i*_nQ+j;

                m_Hp(k, k) += _c/std::sqrt(_w)*xm(i)*xc(j);
            }
        }
    }

    molecular_system<T>& molecule(){return m_molsystem;}
    const molecular_system<T>& molecule() const{return m_molsystem;}

    bool& diagonalise_position() {return m_diagonalise_position;}
    const bool& diagonalise_position() const {return m_diagonalise_position;}

public:
    T mode_frequency() const{return _w;}
protected:  
    linalg::matrix<T> m_Hp;
    linalg::matrix<T> m_R;
    linalg::matrix<T> m_U;
    linalg::vector<T> m_E;
    linalg::vector<T> m_Rv;
    linalg::matrix<T> m_Hpf;
    linalg::matrix<T> temp;
    hermite_dvr<T> m_herm;
    linalg::vector<T> psi0;
    linalg::diagonal_matrix<T> diag;
    molecular_system<T> m_molsystem;
    linalg::eigensolver<linalg::hermitian_matrix<T> > _solver;

    bool m_diagonalise_position;
    size_t _nT;
    size_t _nTmin;
    size_t _nQ;
    T _Ep;

    T _w;
    T _c;

    template <typename Otype>
    void transform_operator(const Otype& op,  linalg::matrix<T>& mop)
    {
        mop.resize(_nT, _nT);
        temp = m_U*op;
        mop = temp*linalg::trans(m_U);
    }
    
public:
    void load(const rapidjson::Value& obj)
    {
        try
        {
            //T from_eV = 1.0/27.2114;
            T from_cmn1 = 1.0/219474.63;
            ASSERT(obj.HasMember("omega"), "Required parameters are not present.");
            ASSERT(obj["omega"].IsNumber(), "Required parameters are not correctly specified.");
            _w = obj["omega"].GetDouble()*from_cmn1;

            ASSERT(obj.HasMember("c"), "Required parameters are not present.");
            ASSERT(obj["c"].IsNumber(), "Required parameters are not correctly specified.");
            _c = obj["c"].GetDouble();

            ASSERT(obj.HasMember("ep"), "Required parameters are not present.");
            ASSERT(obj["ep"].IsNumber(), "Required parameters are not correctly specified.");
            _Ep = obj["ep"].GetDouble();

            ASSERT(obj.HasMember("nt"), "Required parameters are not present.");
            ASSERT(obj["nt"].IsInt(), "Required parameters are not correctly specified.");
            _nTmin = obj["nt"].GetInt64();

            ASSERT(obj.HasMember("nq"), "Required parameters are not present.");
            ASSERT(obj["nq"].IsInt(), "Required parameters are not correctly specified.");
            _nQ = obj["nq"].GetInt64();

            CALL_AND_HANDLE(m_molsystem.load(obj), "Failed to load the molecule system.");
            m_molsystem.use_dvr() = true;

            if(obj.HasMember("diagpos"))
            {
                ASSERT(obj["diagpos"].IsBool(), "Required parameters are not correctly specified.");
                m_diagonalise_position = obj["diagpos"].GetBool();
            }
            else
            {
                m_diagonalise_position = false;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load exponential cutoff from file.");
        } 
    }
};

#endif

