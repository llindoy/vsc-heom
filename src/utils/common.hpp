#ifndef EOS_CONFIG_COMMON_HPP
#define EOS_CONFIG_COMMON_HPP

#include <linalg/utils/linalg_utils.hpp>
#include <linalg/utils/exception_handling.hpp>
#include <memory>

template <typename T, typename ... Args>
std::unique_ptr<T> make_unique(Args&&... args){return std::unique_ptr<T>(new T(std::forward<Args>(args)...));}

template <typename T>
std::shared_ptr<T> make_shared(){return std::shared_ptr<T>(new T());}

template <typename T, typename ... Args>
std::shared_ptr<T> make_shared(Args&&... args)
{
    return std::shared_ptr<T>(new T(std::forward<Args>(args)...));
}


#endif

