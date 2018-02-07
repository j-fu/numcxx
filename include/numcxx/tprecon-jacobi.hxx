///
/// \file tprecon-jacobi.hxx
///
/// Header for  Jacobi preconditioner
/// 
#ifndef TPRECON_JACOBI_HXX
#define TPRECON_JACOBI_HXX

#include "tsparsematrix.hxx"

namespace numcxx
{

  /// Jacobi preconditioner class
  template<typename T> 
  class TPreconJacobi: public  TLinSolver<T>
  {

  public:
    std::shared_ptr< TArray1<T> > pInvDiag;

    /// Create Preconditioner
    TPreconJacobi(const std::shared_ptr<TSparseMatrix<T>> pA);
    TPreconJacobi(TSparseMatrix<T> &A);

    ~TPreconJacobi(){};

    /// Create preconditioner
    static std::shared_ptr<TPreconJacobi<T>> create(const std::shared_ptr<TSparseMatrix<T>> pA);

    /// Perform actual computation preconditioner
    void update(TSparseMatrix<T> &A);

    /// Solve preconditioning system
    void solve( TArray<T> & Sol,  const TArray<T> & Rhs) const;
  };
}

#include "tprecon-jacobi.ixx"

#endif

