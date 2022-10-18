#ifndef HEOM_BATH_CORRELATION_CONTAINER_HPP
#define HEOM_BATH_CORRELATION_CONTAINER_HPP

#include <linalg/linalg.hpp>
#include "bath_correlation_terms.hpp"
#include "heom_truncation.hpp"


//an object for handling the generalised heom expansion.  Here we choose the sigma vector such that 
//\sigma_k = \sqrt(|S_k|) = \sqrt(|\sum_l s_{kl}\phi_l(0)).  This ensures that the terms that increase 
//and decrease the number of excitations are the same as those in the Shi renormalisation for exponential 

template <typename T>
class heom_bath
{
public:
    heom_bath(){}
    heom_bath(const heom_bath& o) = default;
    heom_bath(heom_bath&& o) = default;

    heom_bath& operator=(const heom_bath& o) = default;
    heom_bath& operator=(heom_bath&& o) = default;

    void add_term(const correlation_terms<T>& term)
    {   
        m_terms.push_back(term);
    }
    void add_term(correlation_terms<T>&& term)
    {   
        m_terms.push_back(std::forward<correlation_terms<T>>(term));
    }

    void remove_term(size_t i)
    {
        ASSERT(i < m_terms.size(), "Index out of bounds."); 
        m_terms.erase(m_terms.begin() + i);
    }

    const correlation_terms<T>& operator[](size_t i) const{ASSERT(i < m_terms.size(), "Index out of bounds."); return m_terms[i];}
    correlation_terms<T>& operator[](size_t i){ASSERT(i < m_terms.size(), "Index out of bounds."); return m_terms[i];}
    const correlation_terms<T>& term(size_t i) const{ASSERT(i < m_terms.size(), "Index out of bounds."); return m_terms[i];}
    correlation_terms<T>& term(size_t i){ASSERT(i < m_terms.size(), "Index out of bounds."); return m_terms[i];}

    template <typename trunc_type> 
    typename std::enable_if<std::is_base_of<heom_truncation<T>, trunc_type>::value, void>::type 
    bind_truncation(const trunc_type& trunc)
    {
        m_trunc = std::make_shared<trunc_type>(trunc);
    }

    template <typename trunc_type> 
    typename std::enable_if<std::is_base_of<heom_truncation<T>, trunc_type>::value, void>::type 
    bind_truncation(trunc_type&& trunc)
    {
        m_trunc = std::make_shared<trunc_type>(std::forward<trunc_type>(trunc));
    }

    bool has_truncation() const
    {
        return (m_trunc != nullptr);
    }

    std::shared_ptr<heom_truncation<T>> truncation() const{return m_trunc;}
    bool has_constant_truncation() const
    {
        if(!has_truncation()){return false;}
        return m_trunc->constant_truncation();
    }

    size_t size() const{return m_terms.size();}
protected:

    std::shared_ptr<heom_truncation<T>> m_trunc;
    std::vector<correlation_terms<T>> m_terms;
};

#endif

