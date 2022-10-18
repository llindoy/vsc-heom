#ifndef EOS_CONFIG_IO_HPP
#define EOS_CONFIG_IO_HPP

#include <regex>
#include <string>
#include <algorithm>
#include <iostream>
#include <linalg/linalg.hpp>

#include "common.hpp"
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>


static std::string& remove_whitespace_and_to_lower(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c){return std::tolower(c);});
    str.erase(std::remove_if(str.begin(), str.end(), [](unsigned char c){return (c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '\v' || c == '\f');}), str.end());
    return str;
}

template <typename U, typename V>
void dfs_replace_space_and_capitals_in_key(U& value, V& allocator)
{
    if(value.IsObject())
    {
        for(rapidjson::Value::MemberIterator itr = value.MemberBegin(); itr != value.MemberEnd(); ++itr)
        {
            std::string s(itr->name.GetString());
            remove_whitespace_and_to_lower(s);
            itr->name.SetString(s.c_str(), allocator);
            dfs_replace_space_and_capitals_in_key(itr->value, allocator);
        }
    }
    else if(value.IsArray())
    {
        for(rapidjson::Value::ValueIterator itr = value.Begin(); itr != value.End(); ++itr)
        {
            dfs_replace_space_and_capitals_in_key(*itr, allocator);
        }
    }
}

template <typename T> 
struct complex_from_string
{
public:
    using complex_type = linalg::complex<T>;
    static bool is_valid(const std::string& str)
    {
        std::regex x_regex("^(-)?([0-9]+([.][0-9]*)?|[.][0-9]+)$");   //for matching a purely real number
        std::regex iy_regex("^(-)?([0-9]+([.][0-9]*)?|[.][0-9]+)[ij]$");  //for matching a purely imaginary number of the form yi
        std::regex x_yi_regex("^((-)?([0-9]+([.][0-9]*)?|[.][0-9]+))(([-+])?([0-9]+([.][0-9]*)?|[.][0-9]+)[ij])$");   //for 
        std::regex yi_x_regex("^((-)?([0-9]+([.][0-9]*)?|[.][0-9]+)[ij])(([-+])?([0-9]+([.][0-9]*)?|[.][0-9]+))$");   //for 
        std::smatch match;
        return  (std::regex_match(str, match, x_regex) || std::regex_match(str, match, iy_regex) || std::regex_match(str, match, x_yi_regex) || std::regex_match(str, match, yi_x_regex));
    }
    static complex_type get(const std::string& str)
    {
        T x = 0, y = 0;
        std::regex x_regex("^(-)?([0-9]+([.][0-9]*)?|[.][0-9]+)$");   //for matching a purely real number
        std::regex iy_regex("^(-)?([0-9]+([.][0-9]*)?|[.][0-9]+)[ij]$");  //for matching a purely imaginary number of the form yi
        std::regex x_yi_regex("^((-)?([0-9]+([.][0-9]*)?|[.][0-9]+))(([-+])?([0-9]+([.][0-9]*)?|[.][0-9]+)[ij])$");   //for 
        std::regex yi_x_regex("^((-)?([0-9]+([.][0-9]*)?|[.][0-9]+)[ij])(([-+])?([0-9]+([.][0-9]*)?|[.][0-9]+))$");   //for 
        std::smatch match;
        if(std::regex_match(str, match, x_regex)){x = std::stod(match[0].str());}
        else if(std::regex_match(str, match, iy_regex)){y = (match[1].matched ? -1.0 : 1.0)*std::stod(match[2].str());}
        else if(std::regex_match(str, match, x_yi_regex))
        {
            x = std::stod(match[1].str());
            y = (match[6].matched && match[6].str() == "-" ? -1.0 : 1.0) * std::stod(match[7].str());
        }
        else if(std::regex_match(str, match, yi_x_regex))
        {
            x = (match[6].matched && match[6].str() == "-" ? -1.0 : 1.0) * std::stod(match[7].str());
            y = (match[2].matched && match[2].str() == "-" ? -1.0 : 1.0) * std::stod(match[3].str());
        }
        else{RAISE_EXCEPTION("Unable to obtaind complex value from string.  Invalid string pattern.");}
        
        return complex_type(x, y);
    }
};  //complex_from_string

template <typename T>
struct parse_complex
{
public:
    static_assert(std::is_floating_point<T>::value, "Require a floating point real type when parsing complex numbers.");
    using real_type = T;
    using complex_type = linalg::complex<T>;
public:
    static bool is_complex(const rapidjson::Value& val)
    {
        if(val.IsNumber()){return true;}
        else if(val.IsObject())
        {
            if(val.HasMember("real") && val.HasMember("imag"))
            {
                if(val["real"].IsNumber() && val["imag"].IsNumber()){return true;}
            }
        }
        else if(val.IsString())
        {
            std::string s(val.GetString());
            remove_whitespace_and_to_lower(s);
            return complex_from_string<T>::is_valid(s);
        }
        return false;
    }

    static complex_type get_complex(const rapidjson::Value& val)
    {
        try
        {
            complex_type z;
            if(val.IsNumber())
            {
                real_type x = val.GetDouble();
                z = complex_type(x, 0);
            }
            else if(val.IsObject())
            {
                if(val.HasMember("real") && val.HasMember("imag"))
                {
                    if(val["real"].IsNumber() && val["imag"].IsNumber())
                    {
                        real_type x = val["real"].GetDouble();
                        real_type y = val["imag"].GetDouble();
                        z = complex_type(x, y);
                    }
                    else{RAISE_EXCEPTION("The real and imaginary parts of the object are not both numbers.");}
                }
                else{RAISE_EXCEPTION("The object did not contain both real and imaginary components.");}
            }
            else if(val.IsString())
            {
                std::string s(val.GetString());
                remove_whitespace_and_to_lower(s);
                z = complex_from_string<T>::get(s);
            }
            return z;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to get complex number from rapidjson value object.");
        }
    }

};  //parse_complex

//some classes to help with loading rapid json objects
template <typename T, typename enabled = void> 
class rapidjson_loader;

template <typename T>
class rapidjson_loader<linalg::complex<T>, typename std::enable_if<std::is_floating_point<T>::value, void>::type >
{
public:
    using complex_type = linalg::complex<T>;
public:
    static bool validate(const rapidjson::Value& val){return parse_complex<T>::is_complex(val);}
    static complex_type load(const rapidjson::Value& val){CALL_AND_RETHROW(return parse_complex<T>::get_complex(val));}
    static void load(const rapidjson::Value& val, complex_type& v){CALL_AND_RETHROW(v = parse_complex<T>::get_complex(val));}
};

template <typename T>
class rapidjson_loader<T, typename std::enable_if<std::is_floating_point<T>::value, void>::type >
{
public:
    using real_type = T;
public:
    static bool validate(const rapidjson::Value& val){return val.IsNumber();}
    static real_type load(const rapidjson::Value& val){ASSERT(val.IsNumber(), "Failed to load real_type the input val object is not a number."); return val.GetDouble();}
    static void load(const rapidjson::Value& val, real_type& v){ASSERT(val.IsNumber(), "Failed to load real_type the input val object is not a number."); v = val.GetDouble();}
};


template <typename T>
class rapidjson_loader<T, typename std::enable_if<std::is_integral<T>::value, void>::type >
{
public:
    using int_type = T;
public:
    static bool validate(const rapidjson::Value& val){return val.IsInt64();}
    static int_type load(const rapidjson::Value& val){ASSERT(val.IsInt64(), "Failed to load int_type the input val object is not a number."); return val.GetInt64();}
    static void load(const rapidjson::Value& val, int_type& v){ASSERT(val.IsInt64(), "Failed to load int_type the input val object is not a number."); v = val.GetInt64();}
};


template <>
class rapidjson_loader<bool, void>
{
public:
public:
    static bool validate(const rapidjson::Value& val){return val.IsBool();}
    static bool load(const rapidjson::Value& val){ASSERT(val.IsBool(), "Failed to load bool the input val object is not a number."); return val.GetBool();}
    static void load(const rapidjson::Value& val, bool& v){ASSERT(val.IsBool(), "Failed to load bool the input val object is not a number."); v = val.GetBool();}
};

template <typename T>
class rapidjson_loader<linalg::vector<T>>
{
protected:
    using load_elements = rapidjson_loader<T>;
public:
    static bool validate(const rapidjson::Value& val)
    {
        if(val.IsArray())
        {
            for(size_t i=0; i < val.Size(); ++i)
            {
                if(!load_elements::validate(val[i])){return false;}
            }
            return true;
        }
        else{return load_elements::validate(val);}
    }

    static linalg::vector<T> load(const rapidjson::Value& val)
    {
        linalg::vector<T> ret;
        CALL_AND_RETHROW(load(val, ret));
        return ret;
    }    


    static void load(const rapidjson::Value& val, linalg::vector<T>& ret)
    {
        try
        {
            ASSERT(validate(val), "Invalid rapidjson object.");
            if(val.IsArray())
            {
                ret.resize(val.Size());
                for(size_t i=0; i < ret.size(); ++i)
                {
                    CALL_AND_HANDLE(load_elements::load(val[i], ret(i)), "Failed to load element of vector."); 
                }
            }
            else
            {
                ret.resize(1);
                CALL_AND_HANDLE(load_elements::load(val, ret(0)), "Failed to load element of vector."); 
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load vector from rapidjson value.");
        }
    }
};

template <typename T>
class rapidjson_loader<linalg::diagonal_matrix<T>>
{
protected:
    using load_elements = rapidjson_loader<T>;
public:
    static bool validate(const rapidjson::Value& val)
    {
        if(val.IsArray())
        {
            for(size_t i=0; i < val.Size(); ++i)
            {
                if(!load_elements::validate(val[i])){return false;}
            }
            return true;
        }
        else{return false;}
    }

    static linalg::diagonal_matrix<T> load(const rapidjson::Value& val)
    {
        linalg::diagonal_matrix<T> ret;
        CALL_AND_RETHROW(load(val, ret));
        return ret;
    }    

    static void load(const rapidjson::Value& val, linalg::diagonal_matrix<T>& ret)
    {
        try
        {
            ASSERT(validate(val), "Invalid rapidjson object.");
            ret.resize(val.Size());
            for(size_t i=0; i < ret.size(); ++i)
            {
                CALL_AND_HANDLE(load_elements::load(val[i], ret(i, i)), "Failed to load element of vector."); 
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load diagonal matrix from rapidjson value.");
        }
    }
};


template <typename T>
class rapidjson_loader<linalg::matrix<T>>
{
protected:
    using load_elements = rapidjson_loader<T>;

public:
    static bool validate(const rapidjson::Value& val)
    {
        if(val.IsArray())
        {
            size_t n = val.Size();
            size_t m;
            for(size_t i=0; i < n; ++i)
            {
                if(val[i].IsArray())
                {
                    if(i > 0)
                    {
                        if(val[i].Size() != m){return false;}
                    }
                    else{ m  = val[i].Size();}
                    for(size_t j = 0; j < m; ++j)
                    {
                        if(!load_elements::validate(val[i][j])){return false;}
                    }
                }
                else{return false;}
            }
            return true;
        }
        else{return false;}
    }

    static linalg::matrix<T> load(const rapidjson::Value& val)
    {
        linalg::matrix<T> ret;
        CALL_AND_RETHROW(load(val, ret));
        return ret;
    }    

    static void load(const rapidjson::Value& val, linalg::matrix<T>& ret)
    {
        try
        {
            ASSERT(validate(val), "Invalid rapidjson object.");
            
            size_t n = val.Size();
            size_t m = 0;
            if(n != 0){m = val[0].Size();}
            ret.resize(n, m);
            for(size_t i=0; i < n; ++i)
            {
                for(size_t j=0; j < m; ++j)
                {
                    CALL_AND_HANDLE(load_elements::load(val[i][j], ret(i, j)), "Failed to load element of matrix."); 
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load matrix from rapidjson value.");
        }
    }
};

template <typename T>
class rapidjson_loader<linalg::csr_matrix<T>>
{
protected:
    using load_elements = rapidjson_loader<T>;
    using coo_type = typename linalg::csr_matrix<T>::coo_type;
    using index_type = typename linalg::csr_matrix<T>::index_type;
public:
    static bool validate(const rapidjson::Value& val)
    {
        if(val.IsArray())
        {
            for(size_t i=0; i < val.Size(); ++i)
            {
                if(val[i].IsArray())
                {
                    if(val[i].Size() == 3)
                    {
                        if(!val[i][0].IsInt() || !val[i][1].IsInt() || !load_elements::validate(val[i][2])){return false;}
                    }
                    else{return false;}
                }
                else{return false;}
            }
            return true;
        }
        else if(val.IsObject())
        {
            if(val.HasMember("rowptr") && val.HasMember("colind") && val.HasMember("data")) 
            {
                if(val["rowptr"].IsArray()){for(size_t i = 0; i < val["rowptr"].Size(); ++i){if(!val["rowptr"][i].IsInt()){return false;}}}
                if(val["colind"].IsArray()){for(size_t i = 0; i < val["colind"].Size(); ++i){if(!val["colind"][i].IsInt()){return false;}}}
                if(val["buffer"].IsArray()){for(size_t i = 0; i < val["buffer"].Size(); ++i){if(!load_elements::validate(val["buffer"][i])){return false;}}}
                if(val.HasMember("ncols"))
                {
                    if(!val["ncols"].IsInt()){return false;}
                }
    
                if(val["colind"].Size() != val["buffer"].Size()){return false;}
            }
            else{return false;}
        }
        else{return false;}
        return false;
    }

    static linalg::vector<T> load(const rapidjson::Value& val)
    {
        linalg::csr_matrix<T> ret;
        CALL_AND_RETHROW(load(val, ret));
        return ret;
    }    

    static void load(const rapidjson::Value& val, linalg::csr_matrix<T>& ret, size_t size)
    {
        try
        {
            ASSERT(validate(val), "Invalid rapidjson object.");
            if(val.IsArray())
            {
                index_type max_row = 0;
                index_type max_col = 0;
                coo_type coo(val.Size());
                for(size_t i = 0; i < val.Size(); ++i)
                {
                    index_type row, col;
                    T data;
                    row = val[i][0].GetInt64();
                    col = val[i][1].GetInt64();
                    ASSERT(row >= 0 && col >= 0, "Row and Column Indices must be positive.");
                    CALL_AND_HANDLE(load_elements::load(val[i][2], data), "Failed to load data.");
                    coo[i] = std::make_tuple(row, col, data);
                    if(row > max_row){max_row = row;}
                    if(col > max_col){max_col = col;}
                }

                ASSERT(static_cast<size_t>(max_row) < size && static_cast<size_t>(max_col) < size , "An element is out of bounds");
                ret.init(coo, size, size);
            }
            else if(val.IsObject())
            {
                size_t nrows = val["rowptr"].Size()-1;
                int64_t _nnz = val["rowptr"][nrows].GetInt64();

                ASSERT(nrows < size, "Index out of bound.");
                size_t ncols = size;

                ASSERT(val["buffer"].Size() == _nnz, "Invalid sized buffer.");
                size_t nnz = _nnz;

                ret.resize(nnz, size, size);
    
                auto* buffer = ret.buffer();
                auto* colind = ret.colind();
                auto* rowptr = ret.rowptr();

                //read in the rowptr objects that have been specified
                for(size_t i = 0; i < nrows+1; ++i)
                {
                    index_type rp = val["rowptr"][i].GetInt64();
                    ASSERT(rp >= 0, "Failed to read in csr matrix a rowptr is negative.");
                    rowptr[i] = rp;
                }
                //and pad with the remaining values
                for(size_t i = nrows+1; i < size+1; ++i){rowptr[i] = nnz;}

                for(size_t i = 0; i < nnz; ++i)
                {
                    index_type col = val["colind"][i].GetInt64();
                    ASSERT(col < static_cast<index_type>(ncols) && col >= 0, "Failed to read in csr matrix a column is out of bounds.");
                    colind[i] = col;

                    T d;
                    CALL_AND_HANDLE(load_elements::load(val["buffer"][i], d), "Failed to load buffer.");
                    buffer[i] = d;
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to load csr matrix from rapidjson value.");
        }
    }
};



#endif /*EOS_CONFIG_IO_HPP*/
