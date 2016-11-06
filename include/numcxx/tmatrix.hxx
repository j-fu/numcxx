#ifndef NUMCXX_TMATRIX_H
#define NUMCXX_TMATRIX_H

#include <cmath>
#include <cstdlib>

#include "tarray1.hxx"
#include "tarray2.hxx"

namespace numcxx
{

    template<typename T> 
    class TMatrix: public TArray2<T>, TLinOperator<T>
    {
    public:
        using TArray2<T>::size;
        using TArray2<T>::shape;
        using TArray<T>::operator[];

        TMatrix();
        TMatrix(index n0, index n1);
        TMatrix(index n0, index n1, T*data, std::function<void(T*p)> deleter);
        TMatrix(index n0,index n1, T*data, std::shared_ptr<void> datamanager);
        TMatrix(const  std::initializer_list<std::initializer_list<T>> &il );
        static std::shared_ptr<TMatrix <T> > create(index n0,index n1);
        static std::shared_ptr<TMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);
        std::shared_ptr<TMatrix <T> > copy() const;
        std::shared_ptr<TMatrix <T> > clone() const;

        T xentry(const index i, const index j) const {return _data[_idx(i,j)];}


        template <typename VAL>
        TMatrix<T>&  operator=(const VAL  &expr)  { assign(*this,expr); return *this;}

        TMatrix<T>&  operator=(const TMatrix<T> &expr) { assign(*this,expr); return *this;}

        void apply(const TArray1<T> &u, TArray1<T> &v) const;


    private:
        using TArray2<T>::_data;
        using TArray2<T>::_idx;
        using TArray<T>::_check_square;

    };
    
}

#include "tmatrix-imp.hxx"
#endif
