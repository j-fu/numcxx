#ifndef TSOLVER_LAPACK_LU_HXX
#define TSOLVER_LAPACK_LU_HXX

#include "tmatrix.hxx"

namespace numcxx
{
    template<typename T> 
    class TSolverLapackLU: public  TLinSolver<T>
    {
        const std::shared_ptr< TMatrix<T> > a;
        const std::shared_ptr< TMatrix<T> >lu;
        const std::shared_ptr< TArray1<int>> ipiv;
    public:
        TSolverLapackLU(const std::shared_ptr<TMatrix<T>> a);
        void update();
        void solve( TArray1<T> & sol,  const TArray1<T> & rhs);
        static std::shared_ptr<TSolverLapackLU<T>> create(const std::shared_ptr<TMatrix<T>> a);
    };
}
#include "tsolverlapacklu-imp.hxx"
#endif
