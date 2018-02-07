#include <cassert>
namespace numcxx
{
    

  template<typename T> 
  inline TPreconILU<T>::TPreconILU(const std::shared_ptr<TSparseMatrix<T>> pMatrix)
  {
    pInvDiag=std::make_shared<TArray1<T>>(pMatrix->shape(0));
    pDiagIdx=std::make_shared<TArray1<int>>(pMatrix.shape(0));
    update(*pMatrix);
  }

  template<typename T> 
  inline TPreconILU<T>::TPreconILU(TSparseMatrix<T> & A)
  {
    pInvDiag=std::make_shared<TArray1<T>>(A.shape(0));
    pDiagIdx=std::make_shared<TArray1<int>>(A.shape(0));
    update(A);
  }

    
  template<typename T> 
  inline std::shared_ptr<TPreconILU<T>> TPreconILU<T>::create(const std::shared_ptr<TSparseMatrix<T>> pA)
  {
    return std::make_shared<TPreconILU<T>>(pA);
  }
    
  /// Perform actual computation of LU factorization
  template<typename T> 
  inline void TPreconILU<T>::update(TSparseMatrix<T> &M)
  {

    int n=M.shape(0);
    pA=M.pA;
    pIA=M.pIA;
    pJA=M.pJA;

    auto & A=*pA;
    auto & IA=*pIA;
    auto & JA=*pJA;
    auto & XD=*pInvDiag;
    auto & ID=*pDiagIdx;

    for (int i=0; i<n; i++)
    {
      for (int j=IA(i);j<IA(i+1);j++)
        if (JA(j)==i)
        {
          ID(i)=j;
          break;
        }
      XD(i)=A(ID(i));
    }
    
    for (int i=0; i<n; i++)
    {
      XD(i)=1.0/XD(i);
      for(int j=ID(i)+1;j<IA(i+1);j++)
      {
        int k;
        bool found=false;
        for(k=IA(JA(j));k<IA(JA(j)+1);k++)
        {
          if (JA(k)==i)
          {
            found=true;
            break;
          }
        }
        if (found) XD(JA(j))-=A(k)*XD(i)*A(j);
      }
    }

   M.pattern_changed(false);
  }
    
  /// Solve LU factorized system
  template<typename T> 
  inline void TPreconILU<T>::solve(TArray<T> & Sol,  const TArray<T> & Rhs) const
  {
    int n=pIA->size()-1;
    auto & A=*pA;
    auto & IA=*pIA;
    auto & JA=*pJA;
    auto & XD=*pInvDiag;
    auto & ID=*pDiagIdx;

    Sol.resize(Rhs.size());
    for (int i=0; i<n; i++)
    {
      T x=0.0;
      for(int j=IA(i);j<ID(i);j++)
        x+=A(j)*Sol(JA(j));
      Sol(i)=XD(i)*(Rhs(i)-x);
    }

    for (int i=n-1; i>=0; i--)
    {
      T x=0.0;
      for(int j=ID(i)+1;j<IA(i+1);j++)
        x += Sol(JA(j))*A(j);
      Sol(i)-=x*XD(i);
    }
  }

}
