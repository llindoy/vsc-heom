#ifndef FACTOR_TEST_HPP
#define FACTOR_TEST_HPP

#include "factory.hpp"

namespace eos
{
template <typename T>
class F1
{
public:
    F1(){}
    F1(const F1& o) = default;
    virtual ~F1(){}
    virtual std::shared_ptr<F1<T>> clone() const = 0;
    virtual void load(const rapidjson::Value& obj) = 0;
protected:
    void load(const rapidjson::Value& obj, std::string bname)
    {
    }

    F1(const rapidjson::Value& obj, std::string bname){CALL_AND_RETHROW(load(obj, bname));}
};

template <typename T>
class F2 : public F1<T>, public registered_in_factory<F1<T>, F2<T>>
{
public:
    F2(){}
    F2(const rapidjson::Value& obj){CALL_AND_RETHROW(load(obj));}
    F2(const F2& o) = default;
    virtual ~F2(){}
    std::shared_ptr<F1<T>> clone() const final{return std::make_shared<F2<T>>(*this);}
    void load(const rapidjson::Value& obj){}
protected:
};

REGISTER_TEMPLATE_TYPE_INFO_WITH_NAME(F2, "f2", "f2_desc", "options")
template class F2<double>;       


template <typename T>
class F3 : public F1<T>, public registered_in_factory<F1<T>, F3<T>>
{
public:
    F3(){}
    F3(const rapidjson::Value& obj){CALL_AND_RETHROW(load(obj));}
    F3(const F3& o) = default;
    virtual ~F3(){}
    std::shared_ptr<F1<T>> clone() const final{return std::make_shared<F3<T>>(*this);}
    void load(const rapidjson::Value& obj){}
protected:
};

REGISTER_TEMPLATE_TYPE_INFO_WITH_NAME(F3, "f3", "f3_desc", "options")
template class F3<double>;      //explicit instantiation of the          
template class F3<float>;      //explicit instantiation of the          
}


#endif

