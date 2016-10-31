#ifndef NUMCXX_TMATRIX_H
#define NUMCXX_TMATRIX_H

#include <tuple>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "tarray1.hxx"
#include "tarray2.hxx"

namespace numcxx
{
    template<typename T> 
    class TMatrix: public TArray2<T>
    {
    public:
        using TArray2<T>::size;
        using TArray2<T>::shape;

        TMatrix();
        TMatrix(index n0, index n1);
        TMatrix(index n0, index n1, T*data, std::function<void(T*p)> deleter);
        TMatrix(index n0,index n1, T*data, std::shared_ptr<void> datamanager);
        TMatrix(const  std::initializer_list<std::initializer_list<T>> &il );
        static std::shared_ptr<TMatrix <T> > create(index n0,index n1);
        static std::shared_ptr<TMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);
        std::shared_ptr<TMatrix <T> > copy() const;
        std::shared_ptr<TMatrix <T> > clone() const;

        static void lu_decomp(TMatrix<T> &lu, TArray1<int> &ipiv); 
        static void lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,TArray1<T> &sol, const TArray1<T> &rhs);
        void solve(TArray1<T> &sol, const TArray1<T> &rhs) const;
        void apply(const TArray1<T> &u, TArray1<T> &v);

        static void inverse(TMatrix<T> &lu, TArray1<int> &ipiv, TMatrix<T>& inverse);
        void inverse(TMatrix<T>& inverse) const;

        std::tuple<std::shared_ptr<TMatrix<T>>, std::shared_ptr<TArray1<int>>> lu_decomp() const;
        static std::shared_ptr<TArray1<T>> lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,const TArray1<T> &rhs);
        std::shared_ptr<TArray1<T>> solve(const TArray1<T> &rhs) const;
        std::shared_ptr<TArray1<T>> apply(const TArray1<T> &u);

        static std::shared_ptr<TMatrix<T>> inverse(const TMatrix<T> &lu, const TArray1<int> &ipiv);
        std::shared_ptr<TMatrix<T>> inverse() const;
        
        


    private:
        using TArray2<T>::_data;
        using TArray2<T>::_idx;
        using TArray<T>::_check_square;

    };
    

}

#include "tmatrix-def.hxx"
#endif
