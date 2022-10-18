#ifndef MOLECULAR_BASIS_HPP
#define MOLECULAR_BASIS_HPP

#include "potential.hpp"


//sets up the system Hamiltonian for a Caldeira leggett model
template <typename T> 
class molecular_system
{
public:
    molecular_system(){}
    molecular_system(const rapidjson::Value& obj){CALL_AND_RETHROW(load(obj));}
    molecular_system(size_t np, size_t nm, T Rmax, T mass) : _nprim(np), _nM(nm), _Rmax(Rmax), _mass(mass), _use_dvr(false) {}
    ~molecular_system(){}

    template <typename Ptype>
    void bind_potential(Ptype&& p)
    {
        pot = std::make_shared(std::forward<Ptype>(p));
    }

    template <typename Ptype>
    void bind_potential(const Ptype& p)
    {
        pot = std::make_shared(p);
    }

    const T& mass() const{return _mass;}
    T& mass() {return _mass;}

    const T& Rmax() const{return _Rmax;}
    T& Rmax() {return _Rmax;}
    
    const size_t& Nprim() const{return _nprim;}
    size_t& Nprim() {return _nprim;}

    const size_t& Nm() const{return _nM;}
    size_t& Nm() {return _nM;}

    const bool& use_dvr() const{return _use_dvr;}
    bool& use_dvr() {return _use_dvr;}

    const linalg::vector<T>& x() const{ASSERT(_use_dvr, "Need to use dvr to have grid points."); return m_vec;}
    const linalg::vector<T>& E() const{ASSERT(!_use_dvr, "Need to not be using dvr to return energies."); return m_vec;}

    void construct_basis()
    {
        //first set up the renormalised system Hamiltonian
        size_t nR = _nprim;
        T Rmax = _Rmax;
        temp.resize(nR, nR);
        linalg::matrix<T> Ut(nR, nR);  
        linalg::diagonal_matrix<T> Et(nR);  
        m_vec.resize(nR);
        psi0.resize(_nM);   psi0.fill_zeros();
        form_hamiltonian();

        //and diagonalise it to find the relevant fbr states
        CALL_AND_HANDLE(_solver(Hm, Et, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

        bool use_plus = true;
        T ct, st;

        //now truncate the system to only include the nm states we wish to retain
        m_vec.resize(_nM);  m_U.resize(_nM, _nprim);
        for(size_t i = 0; i < _nM; ++i)
        {
            m_vec(i) = Et(i, i);
            std::cerr << m_vec(i)*27.2114 << std::endl;
            for(size_t j = 0; j < _nprim; ++j)
            {
                m_U(i, j) = Ut(j, i);
            }

            linalg::matrix<T> Rij(_nM, _nM);    Rij.fill_zeros();
            //form the r matrix and determine the linear combination of the two lowest energy eigenstates that maximises the R expectation value
            {
                diag.resize(nR, nR);
                for(size_t j = 0; j<nR; ++j)
                {
                    T R = (j)*(2*Rmax)/(nR-1) - Rmax;
                    diag(j, j) = R;                               
                }
                transform_operator(diag, Rij);
            }

            T theta = 0.5*std::atan2(Rij(0, 1) + Rij(1, 0), Rij(1,1)- Rij(0, 0));


            ct = std::cos(theta);
            st = std::sin(theta);
            _R0 = Rij(0, 0)* ct*ct + Rij(1,1)*st*st + (use_plus ? 1.0 : -1.0)*ct*st*(Rij(0,1)+Rij(1,0));
            if(ct*st*(Rij(0, 1) + Rij(1, 0)) >= 0){use_plus = false;}


            //for psi0 as the left most
            if(!use_plus)
            {
                psi0(0) = ct;   psi0(1) = -st;
            }
            else
            {
                psi0(0) = st;   psi0(1) = ct;
            }
        }

        //following this if we are using the dvr then we represent the position operator in the fbr system and obtain the states that diagonalise this operator
        if(_use_dvr)
        {
            linalg::matrix<T> x;
            //first construct the position operator in the primitive dvr system
            for(size_t i =0; i < _nprim; ++i)
            {
                Et(i, i) = (i)*(2*Rmax)/(nR-1) - Rmax;
            }

            //now we transform to the new system
            transform_operator(Et, x);

            _solver(x, m_vec, m_Ufbr_to_dvr);
            Ut = linalg::trans(m_Ufbr_to_dvr)*m_U;
            m_U = Ut;
            linalg::vector<T> psi2 = linalg::trans(m_Ufbr_to_dvr)*psi0;
            psi0 = psi2;
        }
    }

    void position_operator(linalg::matrix<T>& op, T Rdisp=0)
    {
        op.resize(_nM, _nM);
        op.fill_zeros();
        if(_use_dvr)
        {
            for(size_t i = 0; i < _nM; ++i){op(i, i) = m_vec(i) - Rdisp;}
        }
        else
        {
            size_t nR = _nprim;
            T Rmax = _Rmax;
            diag.resize(nR, nR);
            for(size_t i = 0; i<nR; ++i)
            {
                T R = (i)*(2*Rmax)/(nR-1) - Rmax;
                diag(i, i) = R-Rdisp;                               
            }
            transform_operator(diag, op);
        }
    }

    void side_operator(linalg::matrix<T>& op)
    {
        op.resize(_nM, _nM);
        op.fill_zeros();
        if(_use_dvr)
        {
            for(size_t i = 0; i < _nM; ++i){op(i, i) = m_vec(i) < 0 ? 1.0 : 0.0;}
        }
        else
        {
            size_t nR = _nprim;
            T Rmax = _Rmax;
            diag.resize(nR, nR);
            for(size_t i = 0; i<nR; ++i)
            {
                T R = (i)*(2*Rmax)/(nR-1) - Rmax;
                diag(i, i) = R < 0 ? 1.0 : 0.0;                               
            }
            transform_operator(diag, op);
        }
    }

    void R2_operator(linalg::matrix<T>& op)
    {
        op.resize(_nM, _nM);
        op.fill_zeros();
        if(_use_dvr)
        {
            for(size_t i = 0; i < _nM; ++i){op(i, i) = m_vec(i)*m_vec(i);}
        }
        else
        {
            size_t nR = _nprim;
            T Rmax = _Rmax;
            diag.resize(nR, nR);
            for(size_t i = 0; i<nR; ++i)
            {
                T R = (i)*(2*Rmax)/(nR-1) - Rmax;
                diag(i, i) = R*R;                               
            }
            transform_operator(diag, op);
        }
    }

    void kinetic_energy_operator(linalg::matrix<T>& op)
    {
        //first set up the renormalised system Hamiltonian
        size_t nR = _nprim;
        T Rmax = _Rmax;
        Hm.resize(nR, nR);   Hm.fill_zeros();
        T Ki = M_PI/(2*Rmax/nR);

        //set up the system Hamiltonian in a primitive dvr system
        for(size_t i = 0; i<nR; ++i)
        {
            T R = (i)*(2*Rmax)/(nR-1) - Rmax;
            Hm(i, i) += 0.5/_mass*Ki*Ki/3.0*(1+2.0/(nR*nR));                               //add on the bath renormalisation terms
            for(size_t j=0; j < nR; ++j)
            {
                if(i != j)
                {
                    int64_t diff_ij = static_cast<int64_t>(i) - static_cast<int64_t>(j);
                    Hm(i, j) += 0.5/_mass*2*Ki*Ki/(nR*nR)*std::pow(-1, diff_ij)/(std::sin(diff_ij*M_PI/nR)*std::sin(diff_ij*M_PI/nR));
                }
            }
        }
        transform_operator(Hm, op);
    }

    void hamiltonian_operator(linalg::matrix<T>& op)
    {
        if(_use_dvr)
        {
            form_hamiltonian();
            transform_operator(Hm, op);
        }
        else    
        {
            op.resize(_nM, _nM);
            op.fill_zeros();
            for(size_t i =0; i < _nM; ++i){op(i, i) = m_vec(i);}
        }
    }

    void primitive_grid(linalg::vector<T>& r)
    {
        r.resize(_nprim);
        for(size_t i = 0; i < _nprim; ++i)
        {
            r(i) = (i)*(2*_Rmax)/(_nprim-1) - _Rmax;
        }
    }
    
    const  linalg::vector<T>& Psi0() const{return psi0;}
    const T R0() const{return _R0;}


    void rho0(linalg::matrix<T>& op, T beta, T renorm)
    {
        op.resize(_nM, _nM);
        op.fill_zeros();
        linalg::matrix<T> H, R2;
        hamiltonian_operator(H);
        R2_operator(R2);
        H = H+renorm*R2;

        temp.resize(_nM, _nM);
        linalg::matrix<T> Ut(_nM, _nM);  
        linalg::diagonal_matrix<T> Et(_nM);  

        //and diagonalise it to find the relevant fbr states
        CALL_AND_HANDLE(_solver(H, Et, Ut), "Failed to compute eigenvalues of matter Hamiltonian.");

        Et = elemental_exp(-beta*Et);

        temp = Ut*Et;
        op = temp*linalg::trans(Ut);
    }
protected:
    void form_hamiltonian()
    {
        size_t nR = _nprim;
        T Rmax = _Rmax;
        Hm.resize(nR, nR);   Hm.fill_zeros();
        T Ki = M_PI/(2*Rmax/nR);

        //set up the system Hamiltonian in a primitive dvr system
        for(size_t i = 0; i<nR; ++i)
        {
            T R = (i)*(2*Rmax)/(nR-1) - Rmax;
            T vR = pot->operator()(R);
            Hm(i, i) += 0.5/_mass*Ki*Ki/3.0*(1+2.0/(nR*nR))+ vR;                               
            for(size_t j=0; j < nR; ++j)
            {
                if(i != j)
                {
                    int64_t diff_ij = static_cast<int64_t>(i) - static_cast<int64_t>(j);
                    Hm(i, j) += 0.5/_mass*2*Ki*Ki/(nR*nR)*std::pow(-1, diff_ij)/(std::sin(diff_ij*M_PI/nR)*std::sin(diff_ij*M_PI/nR));
                }
            }
        }
    }

public:
    template <typename Otype>
    void transform_operator(const Otype& op,  linalg::matrix<T>& mop)
    {
        mop.resize(_nM, _nM);
        temp = m_U*op;
        mop = temp*linalg::trans(m_U);
    }

protected:
    linalg::matrix<T> m_U;
    linalg::matrix<T> m_Ufbr_to_dvr;
    linalg::diagonal_matrix<T> diag;
    linalg::matrix<T> Hm;
    linalg::matrix<T> temp;
    linalg::eigensolver<linalg::hermitian_matrix<T> > _solver;
    linalg::vector<T> psi0;
    T _mass;
    std::shared_ptr<potential<T>> pot;
    
    linalg::vector<T> m_vec;
    T _lambda;

    size_t _nprim;
    size_t _nM;
    T _Rmax;
    T _R0;

    bool _use_dvr;


public:
    void load(const rapidjson::Value& obj)
    {
        try
        {
            ASSERT(obj.HasMember("rmax"), "Required parameters are not present.");
            ASSERT(obj["rmax"].IsNumber(), "Required parameters are not correctly specified.");
            _Rmax = obj["rmax"].GetDouble();

            ASSERT(obj.HasMember("mass"), "Required parameters are not present.");
            ASSERT(obj["mass"].IsNumber(), "Required parameters are not correctly specified.");
            _mass = obj["mass"].GetDouble();

            ASSERT(obj.HasMember("nr"), "Required parameters are not present.");
            ASSERT(obj["nr"].IsInt(), "Required parameters are not correctly specified.");
            _nprim = obj["nr"].GetInt64();

            ASSERT(obj.HasMember("nm"), "Required parameters are not present.");
            ASSERT(obj["nm"].IsInt(), "Required parameters are not correctly specified.");
            _nM = obj["nm"].GetInt64();

            if(obj.HasMember("use_dvr"))
            {
                ASSERT(obj["use_dvr"].IsBool(), "Required parameters are not correctly specified.");
                _use_dvr = obj["use_dvr"].GetBool();
            }
            else
            {
                _use_dvr = false;
            }

            ASSERT(obj.HasMember("pot"), "potential not found.");
            ASSERT(obj["pot"].IsObject(), "potential invalid.");
            CALL_AND_HANDLE(pot = factory<potential<T>>::create(obj["pot"]), "Failed to load potential.");
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load exponential cutoff from file.");
        } 
    }
};

#endif

