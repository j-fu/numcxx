#include <cassert>
namespace numcxx
{
    
  /// Create LU factorization class
  template<typename T> 
  inline TPreconJacobi<T>::TPreconJacobi(const std::shared_ptr<TSparseMatrix<T>> pMatrix)
  {
    pInvDiag=std::make_shared<TArray1<T>>(pMatrix->shape(0));
    update(*pMatrix);
  }

  template<typename T> 
  inline TPreconJacobi<T>::TPreconJacobi(TSparseMatrix<T> & A)
  {
    pInvDiag=std::make_shared<TArray1<T>>(A.shape(0));
    update(A);
  }

    
  /// Create LU factorization class
  template<typename T> 
  inline std::shared_ptr<TPreconJacobi<T>> TPreconJacobi<T>::create(const std::shared_ptr<TSparseMatrix<T>> pA)
  {
    return std::make_shared<TPreconJacobi<T>>(pA);
  }
    
  /// Perform actual computation of LU factorization
  template<typename T> 
  inline void TPreconJacobi<T>::update(TSparseMatrix<T> &M)
  {
    int n=M.shape(0);
    assert(n==pInvDiag->size());
    auto &InvDiag=*pInvDiag;
    for (int i=0;i<n;i++)
      InvDiag(i)=1.0/M(i,i);

   M.pattern_changed(false);
  }
    
  /// Solve LU factorized system
  template<typename T> 
  inline void TPreconJacobi<T>::solve( TArray<T> & Sol,  const TArray<T> & Rhs) const
  {
    Sol.resize(Rhs.size());
    int n=pInvDiag->size();
    auto &InvDiag=*pInvDiag;
    for (int i=0;i<n;i++)
      Sol(i)=Rhs(i)*InvDiag(i);
  }


}
