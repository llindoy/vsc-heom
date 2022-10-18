#ifndef EOS_CONFIG_FACTORY_HPP
#define EOS_CONFIG_FACTORY_HPP

#include <memory>
#include <map>
#include <string>
#include <sstream>
#include "io.hpp"
#include "type_info.hpp"


template <typename IB> 
struct factory_info_type
{
    using create_function = std::shared_ptr<IB>(*)();
    using create_function_obj = std::shared_ptr<IB>(*)(const rapidjson::Value&);

    create_function m_f;
    create_function_obj m_fobj;
    std::string m_desc;
    std::string m_req;
    std::string m_name;
};

template <typename INTERFACE_TYPE>
class factory
{
public:
    using info_type = factory_info_type<INTERFACE_TYPE>;

public:
    factory() = delete;
    static bool register_type(const std::string name, const std::string alias, info_type info)
    {
        if(get_info().find(alias) == get_info().end()){get_info()[alias] = info;}   else{return false;}
        if(get_alias().find(name) == get_alias().end()){get_alias()[name] = alias;} else{return false;}
        return true;
    }

    static bool is_registered(const std::string& name){return (get_alias().find(name) != get_alias().end());}

    static bool is_loadable(const rapidjson::Value& obj)  
    {
        try
        {
            ASSERT(obj.IsObject(), "Invalid rapidjson object.");
            ASSERT(obj.HasMember("type"), "Object does not specify a type.");
            ASSERT(obj["type"].IsString(), "The interface type is not a string");
            std::string s(obj["type"].GetString());
            remove_whitespace_and_to_lower(s);

            return (get_alias().find(s) != get_alias().end());
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct object from rapidjson object.");
        }
    }

    static std::shared_ptr<INTERFACE_TYPE> create(const std::string& name)
    {
        std::string s(name);
        remove_whitespace_and_to_lower(s);
        auto alias_it = get_alias().find(s);
        if(alias_it != get_alias().end())
        {
            auto it = get_info().find(alias_it->second);
            if(it != get_info().end())
            {
                return it->second.m_f();
            }   
            else{RAISE_EXCEPTION("Invalid type requested.");}
        }
        else{RAISE_EXCEPTION("Invalid type requested.");}
    }

    template <typename ... Args>
    static std::shared_ptr<INTERFACE_TYPE> create(const rapidjson::Value& obj, Args&& ... args)  
    {
        try
        {
            ASSERT(obj.IsObject(), "Invalid rapidjson object.");
            ASSERT(obj.HasMember("type"), "Object does not specify a type.");
            ASSERT(obj["type"].IsString(), "The interface type is not a string");
            std::string s(obj["type"].GetString());
            remove_whitespace_and_to_lower(s);

            auto alias_it = get_alias().find(s);
            if(alias_it != get_alias().end())
            {
                auto it = get_info().find(alias_it->second);
                if(it != get_info().end())
                {
                    return it->second.m_fobj(obj, std::forward<Args>(args)...);
                }   
                else{RAISE_EXCEPTION("Invalid type requested.");}
            }
            else{RAISE_EXCEPTION("Invalid type requested.");}
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct object from rapidjson object.");
        }
    }

    static std::string description(const std::string& name)
    {
        std::string s(name);
        remove_whitespace_and_to_lower(s);
        auto alias_it = get_alias().find(s);
        if(alias_it != get_alias().end())
        {
            auto it = get_info().find(alias_it->second);
            if(it != get_info().end())
            {
                return it->second.m_desc;
            }   
            else{RAISE_EXCEPTION("Invalid type requested.");}
        }
        else{RAISE_EXCEPTION("Invalid type requested.");}
    }

    static std::string required_inputs(const std::string& name)
    {
        std::string s(name);
        remove_whitespace_and_to_lower(s);
        auto alias_it = get_alias().find(s);
        if(alias_it != get_alias().end())
        {
            auto it = get_info().find(alias_it->second);
            if(it != get_info().end())
            {
                return it->second.m_req;
            }   
            else{RAISE_EXCEPTION("Invalid type requested.");}
        }
        else{RAISE_EXCEPTION("Invalid type requested.");}
    }

    static std::string get_all_info()
    {
        std::ostringstream oss;
        std::string tab("\t");
        std::string doubletab("\t\t");

        std::istringstream sstr;
        std::string line;
        for(auto it = get_info().begin(); it != get_info().end(); ++it)
        {
            oss << it->second.m_name << ": " << std::endl;
            std::string str = it->second.m_desc;
            sstr.clear();  sstr.str(str);
            while(std::getline(sstr, line)){oss << "\t" << line << std::endl;}
            str = it->second.m_req;
            sstr.clear();  sstr.str(str);
            while(std::getline(sstr, line)){oss << "\t\t" << line << std::endl;}
        }
        return oss.str();
    }

private:
    static std::map<std::string, std::string>& get_alias()
    {
        static std::map<std::string, std::string> m_alias;
        return m_alias;
    }

    static std::map<std::string, info_type>& get_info()
    {
        static std::map<std::string, info_type> m_info;
        return m_info;
    }
};


template <typename T> inline void FORCE_INSTANTIATE(T){}

template <typename INTERFACE_TYPE, typename TYPE>
class registered_in_factory
{
public:
    registered_in_factory() {FORCE_INSTANTIATE(m_registered);}
    virtual ~registered_in_factory() {}
private:
    static int m_registered;
};

template <typename INTERFACE_TYPE, class TYPE> 
int registered_in_factory<INTERFACE_TYPE, TYPE>::m_registered = factory<INTERFACE_TYPE>::register_type
(
    type_info<TYPE>::get_name(), 
    type_info<TYPE>::get_alias(), 
    typename factory<INTERFACE_TYPE>::info_type
    {
        type_info<TYPE>::create, 
        type_info<TYPE>::create_from_obj,   
        type_info<TYPE>::get_description(), 
        type_info<TYPE>::get_required_inputs(),
        type_info<TYPE>::get_name()
    }
);


#endif  //EOS_CONFIG_BATH_FACTORY_HPP//

