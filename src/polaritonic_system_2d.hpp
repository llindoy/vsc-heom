#ifndef POLARITONIC_BASIS_2D_HPP
#define POLARITONIC_BASIS_2D_HPP

#include "molecular_system_2D.hpp"

template <typename T> 
class polaritonic_system_2d
{
public:
    polaritonic_system_2d() : m_diagonalise_position(false){}
    polaritonic_system_2d(const rapidjson::Value& obj) : m_diagonalise_position(false){CALL_AND_RETHROW(load(obj));}

    template <typename ... Args> 
    polaritonic_system_2d(size_t np, size_t nf, Args&& ...args) : m_molsystem(std::forward<Args>(args)...), _nP(np), _nF(nf) {}

    const size_t& Npolariton() const{return _nP;}
    size_t& Npolariton(){return _nP;}

    const size_t& NFock() const{return _nF;}
    size_t& NFock(){return _nF;}

    void construct_basis()
    {
        m_molsystem.construct_basis();
        //first set up the renormalised system Hamiltonian
        form_hamiltonian();
        linalg::matrix<T> Ut(m_Hp.shape(0), m_Hp.shape(1));  
        linalg::diagonal_matrix<T> Et(m_Hp.shape(0));  

        //and diagonalise it to find the relevant fbr states
        CALL_AND_HANDLE(_solver(m_Hp, Et, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

        _nP = 0;
        for(size_t i = 0; i < m_Hp.shape(0);++i)
        {
            T cmn1 = 1.0/219474.63;
            if(Et(i, i)/cmn1 < _Ep){++_nP;}
        }
        if(_nP < _nPmin){_nP = _nPmin;}
        m_E.resize(_nP);

        //now truncate the system to only include the nm states we wish to retain
        m_E.resize(_nP);  m_U.resize(_nP, m_Hp.size(0));
        for(size_t i = 0; i < _nP; ++i)
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
            m_R.resize(_nP, _nP);   m_R.fill_zeros();
            size_t nM = m_molsystem.Nm();
            linalg::matrix<T> xm; 
            m_molsystem.R_operator(xm);
            temp2.resize(nM*_nF, nM*_nF);
            for(size_t i = 0; i<nM; ++i)
            {
                for(size_t i2 = 0; i2 < nM; ++i2)
                {
                    for(size_t j = 0; j < _nF; ++j)
                    {
                        temp2(i*_nF+j, i2*_nF+j) = xm(i, i2);                               
                    }
                }
            }
            transform_operator(temp2, m_R);

            //and diagonalise it to find the relevant fbr states
            CALL_AND_HANDLE(_solver(m_R, m_Rv, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

            //now we need to tranform the transformation matrix and the polaritonic Hamiltonian
            {
                linalg::matrix<T> temp3 = linalg::trans(Ut)*m_U;
                m_U = temp3;
            }

            //now set up the Hamiltonian in the new basis
            {
                linalg::diagonal_matrix<T> Ev(m_E.size(),  m_E.size());
                for(size_t i =0; i < m_E.size(); ++i){Ev(i, i) = m_E(i);}
                linalg::matrix<T> temp3 = Ev*Ut;
                m_Hpf = linalg::trans(Ut)*temp3;
            }
        }
    }

    void hamiltonian_operator(linalg::matrix<T>& op, T /* renorm */ = 0)
    {
        op.resize(_nP, _nP);
        if(!m_diagonalise_position)
        {
            op.fill_zeros();
            for(size_t i =0; i < _nP; ++i){op(i, i) = m_E(i);}
        }
        else
        {
            for(size_t i = 0; i < _nP; ++i)
            {
                for(size_t j=0; j < _nP; ++j)
                {
                    op(i, j) = m_Hpf(i, j);
                }
            }
        }
    }

    const linalg::vector<T>& E() const{return m_E;}


    void dipole_operator(linalg::matrix<T>& op)
    {
        if(_coupleR)
        {
            R_operator(op);
        }
        else
        {
            Q_operator(op);
        }
    }

    void R_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        if(!m_diagonalise_position)
        {
            size_t nM = m_molsystem.Nm();
            linalg::matrix<T> xm; 
            m_molsystem.R_operator(xm);
            temp2.resize(nM*_nF, nM*_nF);
            for(size_t i = 0; i<nM; ++i)
            {
                for(size_t i2 = 0; i2 < nM; ++i2)
                {
                    for(size_t j = 0; j < _nF; ++j)
                    {
                        temp2(i*_nF+j, i2*_nF+j) = xm(i, i2);                               
                    }
                }
            }
            transform_operator(temp2, op);
        }
        else
        {
            for(size_t i =0; i < _nP; ++i){op(i, i) = m_Rv(i);}
        }
    }

    void R2_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::matrix<T> xm; 
        m_molsystem.R2_operator(xm);
        temp2.resize(nM*_nF, nM*_nF);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                for(size_t j = 0; j < _nF; ++j)
                {
                    temp2(i*_nF+j, i2*_nF+j) = xm(i, i2);                               
                }
            }
        }
        transform_operator(temp2, op);
    }

    void Q_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::matrix<T> xm; 
        m_molsystem.Q_operator(xm);
        temp2.resize(nM*_nF, nM*_nF);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                for(size_t j = 0; j < _nF; ++j)
                {
                    temp2(i*_nF+j, i2*_nF+j) = xm(i, i2);                               
                }
            }
        }
        transform_operator(temp2, op);
    }

    void Q2_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::matrix<T> xm; 
        m_molsystem.Q2_operator(xm);
        temp2.resize(nM*_nF, nM*_nF);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                for(size_t j = 0; j < _nF; ++j)
                {
                    temp2(i*_nF+j, i2*_nF+j) = xm(i, i2);                               
                }
            }
        }
        transform_operator(temp2, op);
    }

    void qc_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::vector<T> xc = m_herm.x();
        diag.resize(nM*_nF, nM*_nF);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t j = 0; j < _nF; ++j)
            {
                diag(i*_nF+j, i*_nF+j) = xc(j);                               
            }
        }
        transform_operator(diag, op);
    }


    void qc2_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::vector<T> xc = m_herm.x();
        diag.resize(nM*_nF, nM*_nF);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t j = 0; j < _nF; ++j)
            {
                diag(i*_nF+j, i*_nF+j) = xc(j)*xc(j);                               
            }
        }
        transform_operator(diag, op);
    }

    void side_operator(linalg::matrix<T>& op)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();

        size_t nM = m_molsystem.Nm();
        linalg::matrix<T> xm; 
        m_molsystem.side_operator(xm);
        temp2.resize(nM*_nF, nM*_nF);
        for(size_t i = 0; i<nM; ++i)
        {
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                for(size_t j = 0; j < _nF; ++j)
                {
                    temp2(i*_nF+j, i2*_nF+j) = xm(i, i2);                               
                }
            }
        }
        transform_operator(temp2, op);

    }

    void rho0(linalg::matrix<T>& op, T beta, T renorm1, T renorm2, T renorm3)
    {
        op.resize(_nP, _nP);
        op.fill_zeros();
        linalg::matrix<T> H, R2;
        hamiltonian_operator(H);
        R2_operator(R2);
        H = H+renorm1*R2;
        Q2_operator(R2);
        H = H+renorm2*R2;
        qc2_operator(R2);
        H = H+renorm3*R2;

        temp.resize(_nP, _nP);
        linalg::matrix<T> Ut(_nP, _nP);  
        linalg::diagonal_matrix<T> Et(_nP);  

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
        
        m_herm.w0() = _wc;


        m_Hp.resize(nM*_nF, nM*_nF);    m_Hp.fill_zeros();
        linalg::matrix<T> Hm;
        //first form the molecular Hamiltonian
        m_molsystem.hamiltonian_operator(Hm);
        linalg::matrix<T> xM;
        linalg::matrix<T> x2M;
        if(_coupleR)
        {
            m_molsystem.R_operator(xM);
            m_molsystem.R2_operator(x2M);
        } 
        else
        {
            m_molsystem.Q_operator(xM);
            m_molsystem.Q2_operator(x2M);
        }

        //now form the cavity harmonic oscillator Hamiltonian
        m_herm.Ndvr() = _nF;
        m_herm.construct_quadrature();
        linalg::matrix<T> Hc(_nF, _nF);
        m_herm.H0(Hc);
        linalg::vector<T> xc = m_herm.x();

        //here we are working with frequency scaled positions for the cavity mode.  E.g. we are working with
        //x = sqrt(w_c) q.  So make sure we scale all of the coupling terms and self energy terms appropriately

        T chi = _eta*_wc;
        for(size_t i = 0; i < nM; ++i)
        {
            //first add on the purely molecular terms
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                for(size_t j = 0; j < _nF; ++j)
                {
                    //add on the diagonal terms
                    size_t k = i*_nF+j;
                    size_t k2 = i2*_nF+j;
                    m_Hp(k, k2) += Hm(i,  i2);        //the matter terms

                }
            }

            //first add on the purely molecular terms
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                //add on the dipole self energy terms
                for(size_t j = 0; j < _nF; ++j)
                {   
                    size_t k = i*_nF+j;
                    size_t k2 = i2*_nF+j;
                    m_Hp(k, k2) += chi*_eta*x2M(i, i2);        //the dipole self energy terms
                }
            }


            //now we add on the purely cavity terms
            for(size_t j = 0; j < _nF; ++j)
            {
                for(size_t j2 = 0; j2 < _nF; ++j2)
                {
                    //add on the diagonal terms
                    size_t k = i*_nF+j;
                    size_t k2 = i*_nF+j2;
                    m_Hp(k, k2) += Hc(j,  j2);        //the matter terms - including the dipole self energy terms
                }
            }

            //first add on the purely molecular terms
            for(size_t i2 = 0; i2 < nM; ++i2)
            {
                //now add in the molecule cavity coupling
                for(size_t j = 0; j < _nF; ++j)
                {
                    //add on the diagonal terms
                    size_t k = i*_nF+j;
                    size_t k2 = i2*_nF+j;

                    m_Hp(k, k2) += xM(i, i2)*xc(j)*chi*sqrt(2.0);
                }
            }
        }
    }

    molecular_system_2D<T>& molecule(){return m_molsystem;}
    const molecular_system_2D<T>& molecule() const{return m_molsystem;}

public:
    bool& diagonalise_position() {return m_diagonalise_position;}
    const bool& diagonalise_position() const {return m_diagonalise_position;}

public:
    const T& cavity_frequency() const{return _wc;}
    T& cavity_frequency(){return _wc;}

    T mode_frequency() const{return m_molsystem.mode_frequency();}
protected:  
    linalg::matrix<T> m_Hp;
    linalg::matrix<T> m_R;
    linalg::matrix<T> m_U;
    linalg::vector<T> m_E;
    linalg::vector<T> m_Rv;
    linalg::matrix<T> m_Hpf;
    linalg::matrix<T> temp;
    linalg::matrix<T> temp2;
    hermite_dvr<T> m_herm;
    linalg::vector<T> psi0;
    linalg::diagonal_matrix<T> diag;
    molecular_system_2D<T> m_molsystem;
    linalg::eigensolver<linalg::hermitian_matrix<T> > _solver;

    bool m_diagonalise_position;
    size_t _nP;
    size_t _nPmin;
    size_t _nF;
    T _Ep;

    T _wc;
    T _eta;

    bool _coupleR;

    template <typename Otype>
    void transform_operator(const Otype& op,  linalg::matrix<T>& mop)
    {
        mop.resize(_nP, _nP);
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
            ASSERT(obj.HasMember("nu_c"), "Required parameters are not present.");
            ASSERT(obj["nu_c"].IsNumber(), "Required parameters are not correctly specified.");
            _wc = obj["nu_c"].GetDouble()*from_cmn1;

            ASSERT(obj.HasMember("eta"), "Required parameters are not present.");
            ASSERT(obj["eta"].IsNumber(), "Required parameters are not correctly specified.");
            _eta = obj["eta"].GetDouble();

            ASSERT(obj.HasMember("ep"), "Required parameters are not present.");
            ASSERT(obj["ep"].IsNumber(), "Required parameters are not correctly specified.");
            _Ep = obj["ep"].GetDouble();

            ASSERT(obj.HasMember("np"), "Required parameters are not present.");
            ASSERT(obj["np"].IsInt(), "Required parameters are not correctly specified.");
            _nPmin = obj["np"].GetInt64();

            ASSERT(obj.HasMember("nf"), "Required parameters are not present.");
            ASSERT(obj["nf"].IsInt(), "Required parameters are not correctly specified.");
            _nF = obj["nf"].GetInt64();

            CALL_AND_HANDLE(m_molsystem.load(obj), "Failed to load the molecule system.");

            if(obj.HasMember("coupleq"))
            {
                ASSERT(obj["coupleq"].IsBool(), "Required parameters are not correctly specified.");
                _coupleR = !(obj["coupleq"].GetBool());
            }
            else
            {
                _coupleR = true;
            }

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

