#ifndef NUMCXX_TARRAY_H
#define NUMCXX_TARRAY_H

#include <typeinfo>
#include <memory>
#include <stdexcept> 

namespace numcxx
{
    using index= unsigned int;
    
    template<typename T> class TArray
    {
    public:
        T*data() const;
        index ndim() const;
        index size() const;
        index shape(const index dim)  const;
        T & operator()(index i0, index i1);
        T & operator()(index i0);
        void fill(T x);
        void operator=(const T a);
        void operator+=(const T a);
        void operator*=(const T a);
        void operator-=(const T a);
        void operator/=(const T a);
        T min() const;
        T max() const;
        T sum() const;
        T norm2() const;
        T norm1() const;
        T normi() const;
        void fill(const TArray<T> &A);
        void fill(std::function< T(const T)> f,const TArray<T> & A);
        static void operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B);
        static void operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C);
        static T dot(const TArray<T>& A,const TArray<T> &B);


    private:
        std::shared_ptr<void>_datamanager =0;
        void _check_bounds(index acc_dim, index acc_ndim, index acc_idx) const;
        const std::function<void(T*p)> _deleter;
        index _shape[3];
        const index _ndim;
        index _size;

    protected:
        T* _data=nullptr;
        void _check_square() const;
        index _idx(index i0) const;
        index _idx(index i0,index i1)  const;
        index _idx(index i0,index i1,index i2)  const;
    

        TArray(index n0);
        TArray(index n0, T*data, std::function<void(T*p)> deleter);
        TArray(index n0, T*data, std::shared_ptr<void> datamanager);

        TArray(index n0, index n1);
        TArray(index n0, index n1, T*data,std::function<void(T*p)> deleter);
        TArray(index n0, index n1, T*data,std::shared_ptr<void> datamanager);

        TArray();
        ~TArray();

    };
}


#include "tarray-def.hxx"
#endif
