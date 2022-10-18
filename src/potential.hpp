#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include "utils/factory.hpp"
#include "utils/io.hpp"

template <typename T> 
class potential
{
public:
    potential(){}
    virtual ~potential(){}
    virtual T operator()(T x) = 0;
    virtual void load(const rapidjson::Value& obj) = 0;

protected:
protected:
    void load(const rapidjson::Value& obj, std::string bname)
    {
        try
        {
            ASSERT(obj.IsObject(), "Invalid rapidjson object.");
            ASSERT(obj.HasMember("type"), "Object does not specify a spectral_density type.");
            ASSERT(obj["type"].IsString(), "The spectral density type is not a string");
            std::string s(obj["type"].GetString());
            remove_whitespace_and_to_lower(s);
            remove_whitespace_and_to_lower(bname);
            ASSERT(s == bname, "The input spectral density type differs from the type that is being created.");
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load abstract_bath object from file.");
        }
    }
    potential(const rapidjson::Value& obj, std::string bname){CALL_AND_RETHROW(load(obj, bname));}
};

template <typename T> 
class quartic : public potential<T>, public registered_in_factory<potential<T>, quartic<T>>
{
public:
    quartic(){}
    quartic(T a1, T a2)  : _a1(a1), _a2(a2) {}
    quartic(const rapidjson::Value& obj) : potential<T>()
    {
        CALL_AND_HANDLE(load(obj), "Failed to construct debye spectral density object from rapidjson value.");
    }
    
    ~quartic(){}

    void load(const rapidjson::Value& obj) final;
    T operator()(T r) final
    {
        return -_a1*r*r + _a2 * r*r*r*r;
    }
protected:
    T _a1;
    T _a2;
};

REGISTER_TEMPLATE_TYPE_INFO_WITH_NAME(quartic, "quartic", "Quartic potential energy", "")

template <typename value_type> 
void quartic<value_type>::load(const rapidjson::Value& obj)
{
    try
    {
        CALL_AND_HANDLE(potential<value_type>::load(obj, type_info<quartic<value_type> >::get_name()), "Failed to load base type variables.");

        value_type from_cmn1 = 1.0/219474.63;
        if(obj.HasMember("wb") && obj.HasMember("eb"))
        {
            ASSERT(obj["wb"].IsNumber(), "Required parameters are not correctly specified.");
            value_type wb = obj["wb"].GetDouble()*from_cmn1;
            ASSERT(wb >= 0, "Invalid cutoff frequency.");

            ASSERT(obj["eb"].IsNumber(), "Required parameters are not correctly specified.");
            value_type eb = obj["eb"].GetDouble()*from_cmn1;
            ASSERT(eb >= 0, "Invalid cutoff frequency.");

            _a1 = wb*wb/2.0;
            _a2 = _a1*_a1/(4*eb);
        }
        else if(obj.HasMember("a1") && obj.HasMember("a2"))
        {
            ASSERT(obj["a1"].IsNumber(), "Required parameters are not correctly specified.");
            _a1 = obj["a1"].GetDouble();

            ASSERT(obj["a2"].IsNumber(), "Required parameters are not correctly specified.");
            _a2 = obj["a2"].GetDouble();
        }
        else{RAISE_EXCEPTION("Input not recognised.");}
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to load exponential cutoff from file.");
    } 
}
template class quartic<double>;

template <typename T> 
class double_well : public potential<T>, public registered_in_factory<potential<T>, double_well<T>>
{
public:
    double_well(): potential<T>(), m_solver(3,3), m_Vm(3,3), m_Em(3,3), m_Xm(3,3){}
    double_well(T w1, T w2, T wb, T R0, T Eb, T Er, T dS) :double_well() {init(w1, w2, wb, R0, Eb, Er, dS);}
    double_well(const rapidjson::Value& obj) : double_well() 
    {
        CALL_AND_HANDLE(load(obj), "Failed to construct debye spectral density object from rapidjson value.");
    }
    
    ~double_well(){}

    void init(T w1, T w2, T wb, T R0, T Eb, T Er, T dS)
    {
        _w1 = w1;
        _w2 = w2;
        _wb = wb;
        _R0 = R0;
        _Eb = Eb;
        _Er = Er;
        _dS = dS;
    }

    void load(const rapidjson::Value& obj) final;
    T operator()(T r) final
    {
        set_vmat(r);
        m_solver(m_Vm, m_Em, m_Xm, true);        
        return m_Em(0, 0);
    }
protected:
    void set_vmat(T x)
    {
        m_Vm.fill_zeros();
        m_Vm(0, 0) = V0(x, -_R0, _w1);
        m_Vm(1, 1) = _Eb - V0(x, 0, _wb) + 100*(1-s(x, -_R0, 1)) + 100*(1-s(x, _R0, -1));
        m_Vm(2, 2) = V0(x, _R0, _w2) + _Er;

        m_Vm(0, 1) = _dS;
        m_Vm(1, 0) = _dS;

        m_Vm(1, 2) = _dS;
        m_Vm(2, 1) = _dS;
    }

    T V0(T r, T r0, T w){return 0.5*w*w*(r-r0)*(r-r0);}
    T s(T r, T r0, T b){return 0.5*(std::tanh(b*(r-r0))+1);}
    


protected:
    T _w1;
    T _w2;
    T _wb;
    T _R0;
    T _Eb;
    T _Er;
    T _dS;

    linalg::eigensolver<linalg::hermitian_matrix<T> > m_solver;    
    linalg::matrix<T> m_Vm;
    linalg::diagonal_matrix<T> m_Em;
    linalg::matrix<T> m_Xm;
};

REGISTER_TEMPLATE_TYPE_INFO_WITH_NAME(double_well, "double_well", "Double well potential energy", "")

template <typename value_type> 
void double_well<value_type>::load(const rapidjson::Value& obj)
{
    try
    {
        CALL_AND_HANDLE(potential<value_type>::load(obj, type_info<double_well<value_type> >::get_name()), "Failed to load base type variables.");

        value_type from_eV = 1.0/27.2114;

        ASSERT(obj.HasMember("w1"), "Required parameters are not present.");
        ASSERT(obj["w1"].IsNumber(), "Required parameters are not correctly specified.");
        value_type w1 = obj["w1"].GetDouble()*from_eV;

        ASSERT(obj.HasMember("w2"), "Required parameters are not present.");
        ASSERT(obj["w2"].IsNumber(), "Required parameters are not correctly specified.");
        value_type w2 = obj["w2"].GetDouble()*from_eV;

        ASSERT(obj.HasMember("wb"), "Required parameters are not present.");
        ASSERT(obj["wb"].IsNumber(), "Required parameters are not correctly specified.");
        value_type wb = obj["wb"].GetDouble()*from_eV;

        ASSERT(obj.HasMember("r0"), "Required parameters are not present.");
        ASSERT(obj["r0"].IsNumber(), "Required parameters are not correctly specified.");
        value_type r0 = obj["r0"].GetDouble();

        ASSERT(obj.HasMember("eb"), "Required parameters are not present.");
        ASSERT(obj["eb"].IsNumber(), "Required parameters are not correctly specified.");
        value_type eb = obj["eb"].GetDouble()*from_eV;

        ASSERT(obj.HasMember("er"), "Required parameters are not present.");
        ASSERT(obj["er"].IsNumber(), "Required parameters are not correctly specified.");
        value_type er = obj["er"].GetDouble()*from_eV;

        ASSERT(obj.HasMember("ds"), "Required parametdss are not present.");
        ASSERT(obj["ds"].IsNumber(), "Required parametdss are not correctly specified.");
        value_type ds = obj["ds"].GetDouble()*from_eV;

        init(w1, w2, wb, r0, eb, er, ds);
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to load exponential cutoff from file.");
    } 
}
template class double_well<double>;


template <typename T> 
class asymmetric_quartic : public potential<T>, public registered_in_factory<potential<T>, asymmetric_quartic<T>>
{
public:
    asymmetric_quartic(){}
    asymmetric_quartic(T a1, T a2, T a3)  : _a1(a1), _a2(a2), _a3(a3) {}
    asymmetric_quartic(const rapidjson::Value& obj) : potential<T>()
    {
        CALL_AND_HANDLE(load(obj), "Failed to construct debye spectral density object from rapidjson value.");
    }
    
    ~asymmetric_quartic(){}

    void load(const rapidjson::Value& obj) final;
    T operator()(T r) final
    {
        return -_a1*r*r + _a2 * r*r*r*r -_a3*r*r*r;
    }
protected:
    T _a1;
    T _a2;
    T _a3;
};

REGISTER_TEMPLATE_TYPE_INFO_WITH_NAME(asymmetric_quartic, "asymmetric_quartic", "Asymmetric quartic potential energy", "")

template <typename value_type> 
void asymmetric_quartic<value_type>::load(const rapidjson::Value& obj)
{
    try
    {
        CALL_AND_HANDLE(potential<value_type>::load(obj, type_info<asymmetric_quartic<value_type> >::get_name()), "Failed to load base type variables.");

        value_type from_cmn1 = 1.0/219474.63;
        if(obj.HasMember("wb") && obj.HasMember("eb"))
        {
            ASSERT(obj["wb"].IsNumber(), "Required parameters are not correctly specified.");
            value_type wb = obj["wb"].GetDouble()*from_cmn1;
            ASSERT(wb >= 0, "Invalid cutoff frequency.");

            ASSERT(obj["eb"].IsNumber(), "Required parameters are not correctly specified.");
            value_type eb = obj["eb"].GetDouble()*from_cmn1;
            ASSERT(eb >= 0, "Invalid cutoff frequency.");

            ASSERT(obj["b"].IsNumber(), "Required parameters are not correctly specified.");
            _a3 = obj["b"].GetDouble()*from_cmn1;

            _a1 = wb*wb/2.0;
            _a2 = _a1*_a1/(4*eb);
        }
        else if(obj.HasMember("a1") && obj.HasMember("a2"))
        {
            ASSERT(obj["a1"].IsNumber(), "Required parameters are not correctly specified.");
            _a1 = obj["a1"].GetDouble();

            ASSERT(obj["a2"].IsNumber(), "Required parameters are not correctly specified.");
            _a2 = obj["a2"].GetDouble();

            ASSERT(obj["b"].IsNumber(), "Required parameters are not correctly specified.");
            _a3 = obj["b"].GetDouble()*from_cmn1;
        }
        else{RAISE_EXCEPTION("Input not recognised.");}
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to load exponential cutoff from file.");
    } 
}
template class asymmetric_quartic<double>;


#endif

