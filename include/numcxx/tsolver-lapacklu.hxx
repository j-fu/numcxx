///
/// \file tsolver-lapacklu.hxx
///
/// Header for numcxx::TSolverLapackLU
/// 
#ifndef TSOLVER_LAPACK_LU_HXX
#define TSOLVER_LAPACK_LU_HXX

#include "tmatrix.hxx"

namespace numcxx
{
  ///
  /// Lapack LU factorization class
  ///
  /// 
  template<typename T> 
  class TSolverLapackLU: public  TLinSolver<T>
  {
    const std::shared_ptr< TMatrix<T> > pMatrix;
    const std::shared_ptr< TMatrix<T> > pLU;
    const std::shared_ptr< TArray1<int>> pIPiv;
  public:
    /// Object constructor, calls update() to obtain
    /// factorization
    TSolverLapackLU(const std::shared_ptr<TMatrix<T>> pMatrix);

    /// Object constructor, calls update(Matriy) to obtain
    /// factorization
    TSolverLapackLU(const TMatrix<T>& Matrix);

    /// Static wrapper around constructor
    static std::shared_ptr<TSolverLapackLU<T>> create(const std::shared_ptr<TMatrix<T>> pMatrix);

    /// Perform computation of LU factorization using actual
    /// state of matrix.
    ///
    /// Uses [``dgetrf``](http://www.netlib.org/lapack/explore-3.1.1-html/dgetrf.f.html)
    /// from the LAPACK library (for T=double)
    void update();

    /// Perform computation of LU factorization using actual
    /// state of matrix.
    void update(const TMatrix<T>& Matrix);

    /// Solve LU factorized system
    ///
    /// Uses [``dgetrs``](http://www.netlib.org/lapack/explore-3.1.1-html/dgetrs.f.html)
    /// from the LAPACK library  (for T=double)
    void solve( TArray<T> & Sol,  const TArray<T> & Rhs) const;

    /// Solve LU factorized system
    void solve( std::shared_ptr<TArray1<T>> Sol,  const std::shared_ptr<TArray1<T>> Rhs) const {solve(*Sol,*Rhs);};

    /// Calculate inverse of matrix A from its LU factors
    /// \todo move constructor
    std::shared_ptr<TMatrix<T>> calculate_inverse();

    /// Default constructor for swig
    TSolverLapackLU() {};
  };
}
#include "tsolver-lapacklu.ixx"
#endif
