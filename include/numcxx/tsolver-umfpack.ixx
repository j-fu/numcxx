///
/// \file tsolver-umfpack.ixx
///
/// Inline method definitions for class numcxx::TSolverUMFPACK
///

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
  inline TSolverUMFPACK<T>::TSolverUMFPACK(const TSparseMatrix<T> &Matrix):
    pMatrix(nullptr),
    Symbolic(nullptr),
    Numeric(nullptr)
  {
    update(Matrix);
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
    if (pMatrix==nullptr)
      throw std::runtime_error("numcxx: TSolverUMFPACK created without smartpointer");
    update(*pMatrix);
  }
  template<typename T> 
  inline void TSolverUMFPACK<T>::update(const TSparseMatrix<T> &Matrix)
  {
    if (Matrix.empty()) return;
    if (Matrix.pattern_changed())
      throw std::runtime_error("numcxx: forgot flush() after sparse pattern changed");
      
    int n=Matrix.shape(0);
    int nia=Matrix.pIA->size();
    int nja=Matrix.pJA->size();


    // double Control[UMFPACK_CONTROL];
    // for (int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0.0;
    // Control[UMFPACK_PIVOT_TOLERANCE]=pivot_tolerance;
    // Control[UMFPACK_DROPTOL]=drop_tolerance;
    // Control[UMFPACK_PRL]=0;
    // double *control=Control;
    // if (pivot_tolerance<0) control=0;
    double *control=nullptr;
        
    int status;
    if (Matrix.pattern_changed() || Symbolic==nullptr)
    {
      if (Symbolic) umfpack_di_free_symbolic(&Symbolic),Symbolic=nullptr;
      if (Numeric) umfpack_di_free_numeric(&Numeric),Numeric=nullptr;
      status=umfpack_di_symbolic (n, n, Matrix.pIA->data(), Matrix.pJA->data(), Matrix.pA->data(), &Symbolic, 0, 0);
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
    status =umfpack_di_numeric (Matrix.pIA->data(), Matrix.pJA->data(), Matrix.pA->data(), Symbolic, &Numeric, control, 0) ;
    if (status>1)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverUMFPACK::update: umfpack_di_numeric error %d\n",status);
      throw std::runtime_error(errormsg);
    }
    /// copy matrix data for use in solve
    pIA=Matrix.pIA;
    pJA=Matrix.pIA;
    pA=Matrix.pA;
  }
    
  /// Solve LU factorized system
  template<typename T> 
  inline void TSolverUMFPACK<T>::solve( TArray<T> & Sol,  const TArray<T> & Rhs)
  {
    Sol.resize(Rhs.size());
    double *control=nullptr;
    int status;
    status=umfpack_di_solve (UMFPACK_At,pIA->data(), pJA->data(),pA->data(), Sol.data(), Rhs.data(), Numeric, control, 0 ) ;
    if (status>1)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverUMFPACK::update: umfpack_di_solve error %d\n",status);
      throw std::runtime_error(errormsg);
    }
  }


}
