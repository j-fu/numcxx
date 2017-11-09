///
/// \file tsolver-umfpack.hxx
///
/// Header for class numcxx::TSolverUMFPACK
///
#ifndef TSOLVER_UMFPACK_HXX
#define TSOLVER_UMFPACK_HXX

#include "tsparsematrix.hxx"

namespace numcxx
{

  /// Bridge class for using umfpack as solver for vmatrix
  ///    
  /// UMFPACK is a GPL licensed sparse  matrix package by Tim Davis of Texas
  /// A&M university,  see [wikipedia](http://en.wikipedia.org/wiki/UMFPACK)
  /// and its [homepage](http://faculty.cse.tamu.edu/davis/suitesparse.html)
  template<typename T> 
  class TSolverUMFPACK: public  TLinSolver<T>
  {
    /// The corresponding matrix
    const std::shared_ptr< TSparseMatrix<T> > pMatrix;
        
    /// Pointer to symbolic factorization data
    void * Symbolic;
        
    /// Pointer to numeric factorization data
    void * Numeric;

    /// Row pointer of matrix
    std::shared_ptr<TArray1<int>> pIA;
        
    /// Column indices
    std::shared_ptr<TArray1<int>> pJA;  
        
    /// Entries
    std::shared_ptr<TArray1 <T> > pA; 

  public:
    TSolverUMFPACK(){};

    /// Create LU factorization class
    TSolverUMFPACK(const std::shared_ptr<TSparseMatrix<T>> pA);

    /// Create LU factorization class
    TSolverUMFPACK(const TSparseMatrix<T>&A);

    TSolverUMFPACK(const TSolverUMFPACK<T>& A)=delete;

    ~TSolverUMFPACK();
    /// Create LU factorization class

    static std::shared_ptr<TSolverUMFPACK<T>> create(const std::shared_ptr<TSparseMatrix<T>> pA);

    /// Perform actual computation of LU factorization
    void update();

    /// Perform actual computation of LU factorization
    void update(const TSparseMatrix<T>&A);


    /// Solve LU factorized system
    void solve( TArray<T> & Sol,  const TArray<T> & Rhs);

    void solve( std::shared_ptr< TArray<T> > Sol,  const std::shared_ptr<TArray<T> > Rhs) {solve(*Sol,*Rhs);};

  };
}

#include "tsolver-umfpack.ixx"

#endif

