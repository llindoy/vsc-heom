#ifndef SUPEROPERATOR_UTILITIES_HPP
#define SUPEROPERATOR_UTILITIES_HPP

#include <memory>
#include <linalg/linalg.hpp>

namespace superoperator
{

template <typename T> 
static inline void csr_to_dense(const linalg::csr_matrix<T>& S, linalg::matrix<T>& D)
{
    CALL_AND_HANDLE(D.resize(S.size(0), S.size(1)), "Failed to resize dense matrix.");   D.fill_zeros();

    auto rowptr = S.rowptr();
    auto colind = S.colind();
    auto buffer = S.buffer();
    for(size_t i = 0; i < S.size(0); ++i)
    {
        for(size_t j = static_cast<size_t>(rowptr[i]); j < static_cast<size_t>(rowptr[i+1]); ++j)
        {
            D(i, colind[j]) = buffer[j];
        }
    }
}

template <typename T, typename RT>
static inline void dense_to_csr(const linalg::matrix<T>& D, linalg::csr_matrix<T>& S, RT prune_tol = -1) 
{
    CALL_AND_HANDLE(S.resize(D.size(0), D.size(1)), "Failed to resize sparse matrix.");
    auto rowptr = S.rowptr();   rowptr[0] = 0;
    size_t nnz = 0;
    for(size_t i=0; i < D.size(0); ++i)
    {
        for(size_t j=0; j < D.size(1); ++j)
        {
            if(prune_tol <= linalg::abs(D(i, j)) )
            {
                ++nnz;
            }
        }
        rowptr[i+1] = nnz;
    }
    CALL_AND_HANDLE(S.resize(nnz), "Failed to set nnz in sparse matrix.");
    
    auto buffer = S.buffer();
    auto colind = S.colind();

    size_t counter = 0;
    for(size_t i=0; i < D.size(0); ++i)
    {
        for(size_t j=0; j < D.size(1); ++j)
        {
            if(prune_tol <= linalg::abs(D(i, j)) )
            {
                buffer[counter] = D(i, j);
                colind[counter] = j;
                ++counter;
            }
        }
    }
}

class superoperator_utilities
{
public:
    template <typename M1, typename U2>
    static inline typename std::enable_if<linalg::is_number<U2>::value, void>::type 
    commutator(const M1& H, linalg::csr_matrix<U2>& L)
    {
        using T = typename linalg::traits<M1>::value_type;
        CALL_AND_HANDLE(construct_superoperator(H, L, [](const T& a){return a;}, [](const T& a){return -a;}), "Failed to construct commutator.");
    }
        
    template <typename M1, typename U2>
    static inline typename std::enable_if<linalg::is_number<U2>::value, void>::type 
    anticommutator(const M1& H, linalg::csr_matrix<U2>& L)
    {
        using T = typename linalg::traits<M1>::value_type;
        CALL_AND_HANDLE(construct_superoperator(H, L, [](const T& a){return a;}, [](const T& a){return a;}), "Failed to construct anticommutator.");
    }


public: 
    template <typename U1>
    static inline void transpose(const linalg::csr_matrix<U1>& H, linalg::csr_matrix<U1>& Ht)
    {
        auto Hbuffer = H.buffer();
        auto Hcolind = H.colind();
        auto Hrowptr = H.rowptr();  

        //now we need to construct the transpose of our csr matrix.  
        using index_type = typename linalg::csr_matrix<U1>::index_type;

        Ht.resize(H.nnz(), H.size(0)+1, H.size(1));
        auto Htbuffer = Ht.buffer();
        auto Htcolind = Ht.colind();
        auto Htrowptr = Ht.rowptr();  
        std::fill(Htrowptr, Htrowptr+H.size(0)+2, 0);

        for(size_t i = 0; i < H.nnz(); ++i)
        {   
            //figure out how many elements per column
            ++Htrowptr[Hcolind[i] + 2];
        }

        //now we go ahead and make a shifted Htrowptr object
        for(size_t i=1; i < H.size(0)+1; ++i){Htrowptr[i+1] += Htrowptr[i];}

        //now we set the transposed column indices and we add in the number of elements in the row we are currently treating at this stage
        for(size_t r = 0; r < H.size(0); ++r)
        {
            for(index_type ri = Hrowptr[r]; ri < Hrowptr[r+1]; ++ri)
            {
                auto index = Htrowptr[Hcolind[ri]+1]++;
                Htbuffer[index] = Hbuffer[ri];
                Htcolind[index] = r;
            }
        }
        //resize the Ht object so that we have the correct number of rows.
        Ht.resize(H.size(0), H.size(1));
    }

protected:
    template <typename U1, typename U2, typename Func1, typename Func2>
    static inline typename std::enable_if<linalg::is_number<U1>::value && linalg::is_number<U2>::value, void>::type 
    construct_superoperator(const linalg::diagonal_matrix<U1>& H, linalg::csr_matrix<U2>& L, Func1&& f1, Func2&& f2)
    {
        ASSERT(H.size(0) == H.size(1), "Unable to form commutator.  Matrix is not square.");
        size_t N = H.size(0);
        size_t nnz = N*N;

        //start by resizing the csr matrix.  For the commutator this is very straightforward
        L.resize(nnz, N*N, N*N);
        auto buffer = L.buffer();
        auto colind = L.colind();
        auto rowptr = L.rowptr();   rowptr[0] = 0;

        size_t counter = 0;
        for(size_t i=0; i < N; ++i)
        {
            for(size_t j = 0; j < N; ++j)
            {
                buffer[counter] = static_cast<U2>(f1(H(i, i)) + f2(H(j, j)));
                colind[counter] = i*N+j;
                ++counter;
                rowptr[i*N+j+1] = counter;
            }
        }
    }
    

    template <typename U1, typename U2, typename Func1, typename Func2>
    static inline typename std::enable_if<linalg::is_number<U1>::value && linalg::is_number<U2>::value, void>::type 
    construct_superoperator(const linalg::csr_matrix<U1>& H, linalg::csr_matrix<U2>& L, Func1&& f1, Func2&& f2)
    {
        ASSERT(H.size(0) == H.size(1), "Unable to form commutator.  Matrix is not square.");
        using index_type = typename linalg::csr_matrix<U1>::index_type;

        index_type N = H.size(0);
        auto Hbuffer = H.buffer();
        auto Hcolind = H.colind();
        auto Hrowptr = H.rowptr();  

        linalg::csr_matrix<U1> Ht;  transpose(H, Ht);
        auto Htbuffer = Ht.buffer();
        auto Htcolind = Ht.colind();
        auto Htrowptr = Ht.rowptr();  

        size_t nnz = 0;
        for(index_type r1 = 0; r1 < N; ++r1)
        {
            for(index_type l1=0; l1 < N; ++l1)
            {   
                index_type ri = Hrowptr[r1];
                index_type li = Htrowptr[l1];
                //this is the outer contribution
                for(; ri < Hrowptr[r1+1] && Hcolind[ri] < r1; ++ri){++nnz;}
                for(; li < Htrowptr[l1+1] && Htcolind[li] < l1; ++li){++nnz;}

                //this is the bit that isn't working at the moment
                if((Hcolind[ri] == r1 && ri < Hrowptr[r1+1]) || (Htcolind[li] == l1 && li < Htrowptr[l1+1]))
                {
                    if(Hcolind[ri] == r1 && ri < Hrowptr[r1+1]){++ri;}
                    if(Htcolind[li] == l1 && li < Htrowptr[l1+1]){++li;}
                    ++nnz;
                }
                for(; li < Htrowptr[l1+1]; ++li){++nnz;}
                for(; ri < Hrowptr[r1+1]; ++ri){++nnz;}
            }
        }


        //start by resizing the csr matrix.  For the commutator this is very straightforward
        L.resize(nnz, N*N, N*N);

        auto Lbuffer = L.buffer();
        auto Lcolind = L.colind();
        auto Lrowptr = L.rowptr();   Lrowptr[0] = 0;
    
        size_t Lcounter = 0;
        for(index_type r1 = 0; r1 < N; ++r1)
        {
            for(index_type l1=0; l1 < N; ++l1)
            {   
                index_type ri = Hrowptr[r1];
                index_type li = Htrowptr[l1];
                //this is the outer contribution
                for(; ri < Hrowptr[r1+1] && Hcolind[ri] < r1; ++ri)
                {
                    Lcolind[Lcounter] = Hcolind[ri]*N + l1;
                    Lbuffer[Lcounter] = f1(Hbuffer[ri]);
                    ++Lcounter;
                }

                for(; li < Htrowptr[l1+1] && Htcolind[li] < l1; ++li)
                {
                    Lcolind[Lcounter] = r1*N + Htcolind[li];
                    Lbuffer[Lcounter] = f2(Htbuffer[li]);
                    ++Lcounter;
                }

                if((Hcolind[ri] == r1 && ri < Hrowptr[r1+1]) || (Htcolind[li] == l1 && li < Htrowptr[l1+1]))
                {
                    Lcolind[Lcounter] = r1*N + l1;
                    Lbuffer[Lcounter] = 0.0;
                    if(ri <  Hrowptr[r1+1] && Hcolind[ri] == r1)
                    {
                        Lbuffer[Lcounter] += f1(Hbuffer[ri]);
                        ++ri;
                    }
                    if(li < Htrowptr[l1+1] && Htcolind[li] == l1)
                    {
                        Lbuffer[Lcounter] += f2(Htbuffer[li]);
                        ++li;
                    }
                    ++Lcounter;
                    
                }

                for(; li < Htrowptr[l1+1]; ++li)
                {
                    Lcolind[Lcounter] = r1*N + Htcolind[li];
                    Lbuffer[Lcounter] = f2(Htbuffer[li]);
                    ++Lcounter;
                }                
                for(; ri < Hrowptr[r1+1]; ++ri)
                {
                    Lcolind[Lcounter] = Hcolind[ri]*N + l1;
                    Lbuffer[Lcounter] = f1(Hbuffer[ri]);
                    ++Lcounter;
                }
                Lrowptr[r1*N+l1+1] = Lcounter;
            }
        }
    }

    template <typename U1, typename U2, typename Func1, typename Func2>
    static inline typename std::enable_if<linalg::is_number<U1>::value && linalg::is_number<U2>::value, void>::type 
    construct_superoperator(const linalg::matrix<U1>& H, linalg::csr_matrix<U2>& L, Func1&& f1, Func2&& f2)
    {
        ASSERT(H.size(0) == H.size(1), "Unable to form commutator.  Matrix is not square.");
        size_t N = H.size(0);
        size_t nnz = N*N*(2*N-1);

        //start by resizing the csr matrix.  For the commutator this is very straightforward
        L.resize(nnz, N*N, N*N);
        auto buffer = L.buffer();
        auto colind = L.colind();
        auto rowptr = L.rowptr();   rowptr[0] = 0;

        size_t counter = 0;
        for(size_t r1=0; r1 < N; ++r1)
        {
            for(size_t r2 = 0; r2 < N; ++r2)
            {
                for(size_t l1 = 0; l1 < r1; ++l1)
                {
                    buffer[counter] = f1(H(r1, l1));
                    colind[counter] = l1*N+r2;
                    ++counter;
                }
                for(size_t l2 = 0; l2 < r2; ++l2)
                {
                    buffer[counter] = f2(H(l2, r2));
                    colind[counter] = r1*N+l2;
                    ++counter;
                }
                buffer[counter] = f1(H(r1, r1)) + f2(H(r2, r2));
                colind[counter] = r1*N+r2;
                ++counter;

                for(size_t l2 = r2+1; l2 < N; ++l2)
                {
                    buffer[counter] = f2(H(l2, r2));
                    colind[counter] = r1*N+l2;
                    ++counter;
                }
                for(size_t l1 = r1+1; l1 < N; ++l1)
                {
                    buffer[counter] = f1(H(r1, l1));
                    colind[counter] = l1*N+r2;
                    ++counter;
                }
            
                rowptr[r1*N+r2+1] = counter;
            }
        }
    }
};




template <typename T> 
class superoperator
{
public:
public:
    superoperator(){}
    superoperator(const superoperator& o) = default;
    superoperator(superoperator&& o) = default;
    superoperator& operator=(const superoperator& o) = default;
    superoperator& operator=(superoperator&& o) = default;
    virtual ~superoperator(){}

    virtual void commutator(linalg::csr_matrix<T>& mat) const = 0;
    virtual void anticommutator(linalg::csr_matrix<T>& mat) const = 0;

    virtual size_t hilb_dim() const = 0;
    virtual size_t liouv_dim() const = 0;
protected:

};




template <typename T>
class dense_operator_superoperator : public superoperator<T> 
{
public:
    dense_operator_superoperator() : superoperator<T>() {}
    dense_operator_superoperator(const linalg::matrix<T>& m) : m_H(m) {}
    dense_operator_superoperator(linalg::matrix<T>&& m) : m_H(std::move(m)) {}
    dense_operator_superoperator(const dense_operator_superoperator& o) = default;
    dense_operator_superoperator(dense_operator_superoperator&& o) = default;

    dense_operator_superoperator& operator=(const linalg::matrix<T>& m){m_H = m;}
    dense_operator_superoperator& operator=(linalg::matrix<T>&& m){ m_H = std::move(m);}
    dense_operator_superoperator& operator=(const dense_operator_superoperator& o) = default;
    dense_operator_superoperator& operator=(dense_operator_superoperator&& o) = default;

    ~dense_operator_superoperator() {}

    virtual void commutator(linalg::csr_matrix<T>& mat) const final
    {   
        superoperator_utilities::commutator(m_H, mat);
    }
    virtual void anticommutator(linalg::csr_matrix<T>& mat) const final
    {
        superoperator_utilities::anticommutator(m_H, mat);
    }

    virtual size_t hilb_dim() const final{ASSERT(m_H.shape(0) == m_H.shape(1), "Invalid matrix shape. Superoperator expects a square matrix."); return m_H.shape(0);}
    virtual size_t liouv_dim() const final{ASSERT(m_H.shape(0) == m_H.shape(1), "Invalid matrix shape. Superoperator expects a square matrix."); return m_H.shape(0)* m_H.shape(0);}
protected:
    linalg::matrix<T> m_H;
};


template <typename T>
class sparse_operator_superoperator : public superoperator<T> 
{
public:
    sparse_operator_superoperator() : superoperator<T>() {}
    sparse_operator_superoperator(const linalg::csr_matrix<T>& m) : m_H(m) {}
    sparse_operator_superoperator(linalg::csr_matrix<T>&& m) : m_H(std::move(m)) {}
    sparse_operator_superoperator(const sparse_operator_superoperator& o) = default;
    sparse_operator_superoperator(sparse_operator_superoperator&& o) = default;

    sparse_operator_superoperator& operator=(const linalg::csr_matrix<T>& m){m_H = m;}
    sparse_operator_superoperator& operator=(linalg::csr_matrix<T>&& m){ m_H = std::move(m);}
    sparse_operator_superoperator& operator=(const sparse_operator_superoperator& o) = default;
    sparse_operator_superoperator& operator=(sparse_operator_superoperator&& o) = default;

    ~sparse_operator_superoperator() {}

    virtual void commutator(linalg::csr_matrix<T>& mat) const final
    {   
        superoperator_utilities::commutator(m_H, mat);
    }
    virtual void anticommutator(linalg::csr_matrix<T>& mat) const final
    {
        superoperator_utilities::anticommutator(m_H, mat);
    }
    virtual size_t hilb_dim() const final{ASSERT(m_H.shape(0) == m_H.shape(1), "Invalid matrix shape. Superoperator expects a square matrix."); return m_H.shape(0);}
    virtual size_t liouv_dim() const final{ASSERT(m_H.shape(0) == m_H.shape(1), "Invalid matrix shape. Superoperator expects a square matrix."); return m_H.shape(0)* m_H.shape(0);}
protected:
    linalg::csr_matrix<T> m_H;
};

template <typename T>
class diagonal_operator_superoperator : public superoperator<T> 
{
    using complex_type = linalg::complex<T>;
public:
    diagonal_operator_superoperator() : superoperator<T>() {}
    diagonal_operator_superoperator(const linalg::diagonal_matrix<T>& m) : m_H(m) {}
    diagonal_operator_superoperator(linalg::diagonal_matrix<T>&& m) : m_H(std::move(m)) {}
    diagonal_operator_superoperator(const diagonal_operator_superoperator& o) = default;
    diagonal_operator_superoperator(diagonal_operator_superoperator&& o) = default;

    diagonal_operator_superoperator& operator=(const linalg::diagonal_matrix<T>& m){m_H = m;}
    diagonal_operator_superoperator& operator=(linalg::diagonal_matrix<T>&& m){ m_H = std::move(m);}
    diagonal_operator_superoperator& operator=(const diagonal_operator_superoperator& o) = default;
    diagonal_operator_superoperator& operator=(diagonal_operator_superoperator&& o) = default;

    ~diagonal_operator_superoperator() {}

    virtual void commutator(linalg::csr_matrix<T>& mat) const final
    {   
        superoperator_utilities::commutator(m_H, mat);
    }
    virtual void anticommutator(linalg::csr_matrix<T>& mat) const final
    {
        superoperator_utilities::anticommutator(m_H, mat);
    }

    void commutator(linalg::diagonal_matrix<T>& mat) const
    {   
        mat.resize(m_H.shape(0)*m_H.shape(0), m_H.shape(1)*m_H.shape(1));   
        for(size_t i =0; i < m_H.shape(0); ++i)
        {
            for(size_t j = 0; j < m_H.shape(1); ++j)
            {
                mat(i*m_H.shape(1)+j, i*m_H.shape(1)+j) = m_H(i,i) - m_H(j, j);
            }
        }
    }
    void anticommutator(linalg::diagonal_matrix<T>& mat) const
    {
        mat.resize(m_H.shape(0)*m_H.shape(0), m_H.shape(1)*m_H.shape(1));   
        for(size_t i =0; i < m_H.shape(0); ++i)
        {
            for(size_t j = 0; j < m_H.shape(1); ++j)
            {
                mat(i*m_H.shape(1)+j, i*m_H.shape(1)+j) = m_H(i,i) + m_H(j, j);
            }
        }
    }

    virtual size_t hilb_dim() const final{ASSERT(m_H.shape(0) == m_H.shape(1), "Invalid matrix shape. Superoperator expects a square matrix."); return m_H.shape(0);}
    virtual size_t liouv_dim() const final{ASSERT(m_H.shape(0) == m_H.shape(1), "Invalid matrix shape. Superoperator expects a square matrix."); return m_H.shape(0)* m_H.shape(0);}
protected:
    linalg::diagonal_matrix<T> m_H;
};


class factory
{
public:
    template <typename T> 
    static std::shared_ptr<superoperator<T>> construct(const linalg::matrix<T>& mat)
    {
        return std::make_shared<dense_operator_superoperator<T> >(mat);
    }

    template <typename T> 
    static std::shared_ptr<superoperator<T>> construct(const linalg::csr_matrix<T>& mat)
    {
        return std::make_shared<sparse_operator_superoperator<T> >(mat);
    }

    template <typename T> 
    static std::shared_ptr<superoperator<T>> construct(const linalg::diagonal_matrix<T>& mat)
    {
        return std::make_shared<diagonal_operator_superoperator<T> >(mat);
    }
};

}   //namespace superoperator

#endif

