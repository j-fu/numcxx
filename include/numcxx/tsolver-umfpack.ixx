#include <cstdio>
#include <umfpack.h>
namespace numcxx
{
    
  /// Create LU factorization class
  template<typename T> 
  inline TSolverUMFPACK<T>::TSolverUMFPACK(const std::shared_ptr<TSparseMatrix<T>> pMatrix):
    pMatrix(pMatrix),
    Symbolic(nullptr),
    Numeric(nullptr)
  {
    update();
  }

  template<typename T> 
  inline TSolverUMFPACK<T>::~TSolverUMFPACK()
  {
    if (Symbolic!=nullptr) 
    {
      umfpack_di_free_symbolic(&Symbolic);
      Symbolic=nullptr;
    }
    if (Numeric!=nullptr)
    {
      umfpack_di_free_numeric(&Numeric);
      Numeric=nullptr;
    }
  }
    
  /// Create LU factorization class
  template<typename T> 
  inline std::shared_ptr<TSolverUMFPACK<T>> TSolverUMFPACK<T>::create(const std::shared_ptr<TSparseMatrix<T>> pA)
  {
    return std::make_shared<TSolverUMFPACK<T>>(pA);
  }
    
  /// Perform actual computation of LU factorization
  template<typename T> 
  inline void TSolverUMFPACK<T>::update()
  {
    pMatrix->flush();
    if (pMatrix->empty()) return;
    int n=pMatrix->shape(0);
    int nia=pMatrix->pIA->size();
    int nja=pMatrix->pJA->size();


    // double Control[UMFPACK_CONTROL];
    // for (int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0.0;
    // Control[UMFPACK_PIVOT_TOLERANCE]=pivot_tolerance;
    // Control[UMFPACK_DROPTOL]=drop_tolerance;
    // Control[UMFPACK_PRL]=0;
    // double *control=Control;
    // if (pivot_tolerance<0) control=0;
    double *control=nullptr;
        
    int status;
    if (pMatrix->pattern_changed() || Symbolic==nullptr)
    {
      if (Symbolic) umfpack_di_free_symbolic(&Symbolic),Symbolic=nullptr;
      if (Numeric) umfpack_di_free_numeric(&Numeric),Numeric=nullptr;
      status=umfpack_di_symbolic (n, n, pMatrix->pIA->data(), pMatrix->pJA->data(), pMatrix->pA->data(), &Symbolic, 0, 0);
      if (status>1)
      {
        char errormsg[80];
        snprintf(errormsg,80,"numcxx::TSolverUMFPACK::update: umfpack_di_symbolic error %d\n",status);
        throw std::runtime_error(errormsg);
      }

    }
    if (Numeric!=nullptr) 
    {
      umfpack_di_free_numeric(&Numeric);
      Numeric=nullptr;
    }
    status =umfpack_di_numeric (pMatrix->pIA->data(), pMatrix->pJA->data(), pMatrix->pA->data(), Symbolic, &Numeric, control, 0) ;
    if (status>1)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverUMFPACK::update: umfpack_di_numeric error %d\n",status);
      throw std::runtime_error(errormsg);
    }
    pMatrix->pattern_changed(false);
  }
    
  /// Solve LU factorized system
  template<typename T> 
  inline void TSolverUMFPACK<T>::solve( TArray<T> & Sol,  const TArray<T> & Rhs)
  {
    Sol.resize(Rhs.size());
    double *control=nullptr;
    int status;
    status=umfpack_di_solve (UMFPACK_At,pMatrix->pIA->data(), pMatrix->pJA->data(), pMatrix->pA->data(), Sol.data(), Rhs.data(), Numeric, control, 0 ) ;
    if (status>1)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverUMFPACK::update: umfpack_di_solve error %d\n",status);
      throw std::runtime_error(errormsg);
    }
  }


}
