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
    template<typename T> class TLinSolver
    {
    public:
        TLinSolver(){};
        virtual void solve( TArray<T> & sol,  const TArray<T> & rhs){};
        virtual void update(void){};
    };

    template<typename T> class TLinOperator
    {
    public:
        TLinOperator(){};
        virtual void apply( const TArray<T> & sol,   TArray<T> & rhs){};
    };



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

        TMatrix<T>&  operator=(const T &expr) { assign(*this,expr); return *this;}

        TMatrix<T>&  operator=(const TMatrix<T> &expr) { assign(*this,expr); return *this;}

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
    
    extern "C"
    {
        extern void dgetrf_(int *n, int *m, double *a, int *lda, int* ipiv, int *info);
        extern void dgetrs_(char *trans,int *n, const int *nrhs, double*a, int* lda, int *ipiv , double *b, int *ldb, int *info );
        extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  double *ALPHA, const double *A, int *LDA,  double *B, int *LDB,  double *BETA,double * C, int *LDC);
        extern void dgetri_(int *n, double*a, int* lda, int *ipiv , double *work, int *lwork, int *info );
    }

    template<typename T> 
    class TSolverLapackLU: public  TLinSolver<T>
    {
        const std::shared_ptr< TMatrix<T> >a;
        const std::shared_ptr< TMatrix<T> >lu;
        const std::shared_ptr< TArray1<int>> ipiv;
    public:
        TSolverLapackLU(const std::shared_ptr<TMatrix<T>> a):
            TLinSolver<T>(),
            a(a),
            lu(a->clone()),
            ipiv(TArray1<int>::create(a->shape(0)))
        {;}
        void update()
        {
            int n=lu->shape(0);
            int info;
            lu->fill(*a);
            dgetrf_(&n,&n,lu->data(),&n,ipiv->data(),&info);
        }
        void solve( TArray<T> & sol,  const TArray<T> & rhs)
        {
            sol.fill(rhs);
            char trans[2]={'T','\0'};
            int n=lu->shape(0);
            int one=1;
            int info;
            dgetrs_(trans,&n,&one,lu->data(),&n,ipiv->data(),sol.data(),&n,&info);
        }

    };


}

#include "tmatrix-imp.hxx"
#endif
