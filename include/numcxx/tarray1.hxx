#ifndef NUMCXX_TARRAY1_H
#define NUMCXX_TARRAY1_H

#include <vector>
#include <ostream>

#include "tarray.hxx"


namespace numcxx
{

    template<typename T>
    class TArray1;


    template<typename T>
    inline std::ostream & operator << (std::ostream & s, TArray1<T> &A);

    template<typename T> class TArray1: public TArray<T>
    {
    public:
        using TArray<T>::size;


        TArray1();
        TArray1(index n0);
        TArray1(index n0, T*data,std::function<void(T*p)> deleter);
        TArray1(index n0, T*data, std::shared_ptr<void> datamanager);
        TArray1(std::shared_ptr<std::vector<T>> v);
        TArray1(const std::initializer_list<T> &il );
        static std::shared_ptr<TArray1 <T> > create(index n1);
        static std::shared_ptr<TArray1 <T> > create(const std::initializer_list<T> il);

        std::shared_ptr<TArray1 <T> > copy() const;
        std::shared_ptr<TArray1 <T> > clone() const;

        T & operator[](index i0);

        T item(index i0);
        void itemset(index i0, T x);

        T __getitem__(index i) const;
        void __setitem__(index i,T d);

        friend std::ostream & operator<< <T>(std::ostream & s, TArray1<T> &A);

        
    private:
        using TArray<T>::_data;
        using TArray<T>::_idx;

    };
}

#include "tarray1-def.hxx"
#endif
