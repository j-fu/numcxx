#ifndef NUMCXX_TARRAY2_H
#define NUMCXX_TARRAY2_H

#include <ostream>
#include "tarray.hxx"


namespace  numcxx
{
    template<typename T>
    class TArray2;

    template<typename T>
    inline std::ostream & operator << (std::ostream & s, TArray2<T> &A);

    template<typename T> 
    class TArray2: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;
        
        TArray2();
        TArray2(index n0, index n1);
        TArray2(index n0, index n1, T*data, std::function<void(T*p)> deleter);
        TArray2(index n0,index n1, T*data, std::shared_ptr<void> datamanager);
        TArray2(const  std::initializer_list<std::initializer_list<T>> &il );
        static std::shared_ptr<TArray2 <T> > create(index n0,index n1);
        static std::shared_ptr<TArray2 <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);
        std::shared_ptr<TArray2 <T> > copy() const;
        std::shared_ptr<TArray2 <T> > clone() const;

        T * operator[](index i0);
        T item(index i0,index i1);
        void itemset(index i0, index i1, T x);
        std::shared_ptr<TArray1 <T> > const __getitem__(index i0);


        friend std::ostream & operator<< <T>(std::ostream & s, TArray2<T> &A);

    protected:
        using TArray<T>::_data;
        using TArray<T>::_idx;
    };

}
#include "tarray2-def.hxx"

#endif
