#ifndef TSOLVER_LAPACK_LU_HXX
#define TSOLVER_LAPACK_LU_HXX

#include "tmatrix.hxx"

namespace numcxx
{
    /// Lapack LU factorization class
    /// 
    /// Th
    template<typename T> 
    class TSolverLapackLU: public  TLinSolver<T>
    {
        const std::shared_ptr< TMatrix<T> > a;
        const std::shared_ptr< TMatrix<T> >lu;
        const std::shared_ptr< TArray1<int>> ipiv;
    public:
        /// Object constructor, calls uptdate() to obtain
        /// factorization
        TSolverLapackLU(const std::shared_ptr<TMatrix<T>> pA);

        /// Static wrapper around constructor
        static std::shared_ptr<TSolverLapackLU<T>> create(const std::shared_ptr<TMatrix<T>> a);

        /// Perform computation of LU factorization using actual
        /// state of matrix.
        ///
        /// Uses [``dgetrf``](http://www.netlib.org/lapack/explore-3.1.1-html/dgetrf.f.html)
        /// from the LAPACK library (for T=double)
        void update();

        /// Solve LU factorized system
        ///
        /// Uses [``dgetrs``](http://www.netlib.org/lapack/explore-3.1.1-html/dgetrs.f.html)
        /// from the LAPACK library  (for T=double)
        void solve( TArray<T> & sol,  const TArray<T> & rhs) const;
    };
}
#include "tsolver-lapacklu.ixx"
#endif
