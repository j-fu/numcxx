///
/// \file tprecon-jacobi.hxx
///
/// Header for  ILU preconditioner
/// 
#ifndef TPRECON_ILU_HXX
#define TPRECON_ILU_HXX

#include "tsparsematrix.hxx"

namespace numcxx
{

  /// ILU preconditioner class
  template<typename T> 
  class TPreconILU: public  TLinSolver<T>
  {

  public:
    std::shared_ptr< TArray1<T> > pInvDiag;
    std::shared_ptr< TArray1<int> > pDiagIdx;

    /// Create Preconditioner
    TPreconILU(const std::shared_ptr<TSparseMatrix<T>> pA);
    TPreconILU(TSparseMatrix<T> &A);

    ~TPreconILU(){};

    /// Create preconditioner
    static std::shared_ptr<TPreconILU<T>> create(const std::shared_ptr<TSparseMatrix<T>> pA);

    /// Perform actual computation preconditioner
    void update(TSparseMatrix<T> &A);

    /// Solve preconditioning system
    void solve(TArray<T> & Sol,  const TArray<T> & Rhs) const;

  private:
    std::shared_ptr<TArray1<T>> pA;
    std::shared_ptr<TArray1<int>> pIA;
    std::shared_ptr<TArray1<int>> pJA;
  };
}

#include "tprecon-ilu.ixx"

#endif

