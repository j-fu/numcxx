#ifndef TSOLVER_LAPACK_LU_HXX
#define TSOLVER_LAPACK_LU_HXX

#include "tmatrix.hxx"

namespace numcxx
{
    /// Lapack LU factorization class
    template<typename T> 
    class TSolverLapackLU: public  TLinSolver<T>
    {
        const std::shared_ptr< TMatrix<T> > a;
        const std::shared_ptr< TMatrix<T> >lu;
        const std::shared_ptr< TArray1<int>> ipiv;
    public:
        /// Create LU factorization class
        TSolverLapackLU(const std::shared_ptr<TMatrix<T>> a);
        /// Create LU factorization class
        static std::shared_ptr<TSolverLapackLU<T>> create(const std::shared_ptr<TMatrix<T>> a);
        /// Perform actual computation of LU factorization
        void update();
        /// Solve LU factorized system
        void solve( TArray1<T> & sol,  const TArray1<T> & rhs);
    };
}
#include "tsolverlapacklu-imp.hxx"
#endif
