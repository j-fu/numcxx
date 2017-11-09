///
/// \file tsolver-lapacklu.ixx
///
/// Inline methods for  numcxx::TSolverLapackLU
/// 
namespace numcxx
{
  template<typename T> 
  inline TSolverLapackLU<T>::TSolverLapackLU(const std::shared_ptr<TMatrix<T>> pMatrix):
    TLinSolver<T>(),
    pMatrix(pMatrix),
    pLU(pMatrix->clone()),
    pIPiv(TArray1<int>::create(pMatrix->shape(0)))
  {
    update();;
  }
  
  template<typename T> 
  inline TSolverLapackLU<T>::TSolverLapackLU(const TMatrix<T> &Matrix):
    TLinSolver<T>(),
    pMatrix(0),
    pLU(TMatrix<T>::create(Matrix.shape(0),Matrix.shape(1))),
    pIPiv(TArray1<int>::create(Matrix.shape(0)))
  { 
    update(Matrix);
  }
    
  template<typename T> 
  inline std::shared_ptr<TSolverLapackLU<T>> TSolverLapackLU<T>::create(const std::shared_ptr<TMatrix<T>> a)
  {
    return std::make_shared<TSolverLapackLU<T>>(a);
  }

  // Declarations for Fortran methos from LAPACK
  extern "C" 
  {
    void sgetrf_(int *n, int *m, float *a, int *lda, int* ipiv, int *info);
    void sgetrs_(char *trans,int *n, const int *nrhs, float*a, int* lda, int *ipiv , float *b, int *ldb, int *info );
    void sgetri_(int *n, float*a, int* lda, int *ipiv , float *work, int *lwork, int *info );
    void dgetrf_(int *n, int *m, double *a, int *lda, int* ipiv, int *info);
    void dgetrs_(char *trans,int *n, const int *nrhs, double*a, int* lda, int *ipiv , double *b, int *ldb, int *info );
    void dgetri_(int *n, double*a, int* lda, int *ipiv , double *work, int *lwork, int *info );
  }

  template<> 
  inline void TSolverLapackLU<double>::update(const TMatrix<double>& Matrix)
  {
    int n=pLU->shape(0);
    int info=0;
    *pLU=Matrix;
    dgetrf_(&n,&n,pLU->data(),&n,pIPiv->data(),&info);
    if (info!=0)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: dgetrf error %d\n",info);
      throw std::runtime_error(errormsg);
    }
  }
    
  template<> 
  inline void TSolverLapackLU<double>::solve( TArray<double> & sol,  const TArray<double> & rhs) const
  {
    assign(sol,rhs);
    char trans[2]={'T','\0'};
    int n=pLU->shape(0);
    int one=1;
    int info=0;
    dgetrs_(trans,&n,&one,pLU->data(),&n,pIPiv->data(),sol.data(),&n,&info);
    if (info!=0)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: dgetrs error %d\n",info);
      throw std::runtime_error(errormsg);
    }
  }


  template<> 
  inline void TSolverLapackLU<float>::update(const TMatrix<float>& Matrix)
  {
    int n=pLU->shape(0);
    int info;
    *pLU=Matrix;
    sgetrf_(&n,&n,pLU->data(),&n,pIPiv->data(),&info);
    if (info!=0)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: sgetrf error %d\n",info);
      throw std::runtime_error(errormsg);
    }

  }

  template<typename T> 
  void TSolverLapackLU<T>::update()
  {
    if (pMatrix==nullptr)
      throw std::runtime_error("numcxx: TSolverLapackLU created without smartpointer");
    update(*pMatrix);
  }
  
    
  template<> 
  inline void TSolverLapackLU<float>::solve( TArray<float> & sol,  const TArray<float> & rhs) const
  {
    assign(sol,rhs);
    char trans[2]={'T','\0'};
    int n=pLU->shape(0);
    int one=1;
    int info;
    sgetrs_(trans,&n,&one,pLU->data(),&n,pIPiv->data(),sol.data(),&n,&info);
    if (info!=0)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: sgetrs error %d\n",info);
      throw std::runtime_error(errormsg);
    }

  }



  template<> 
  inline std::shared_ptr<TMatrix<double>> TSolverLapackLU<double>::calculate_inverse()
  {
    int n=pLU->shape(0);
    auto pInverse=std::make_shared<TMatrix<double>>(*pLU);
    int info;
    TArray1<double> Work(n);
    dgetri_(&n, 
            pInverse->data(), 
            &n, 
            pIPiv->data(),
            Work.data(),
            &n,
            &info);
    if (info!=0)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverLapackLU::calculate_inverse: dgetri error %d\n",info);
      throw std::runtime_error(errormsg);
    }
    return pInverse;
  }

  template<> 
  inline std::shared_ptr<TMatrix<float>> TSolverLapackLU<float>::calculate_inverse()
  {
    int n=pLU->shape(0);
    auto pInverse=std::make_shared<TMatrix<float>>(*pLU);
    int info;
    TArray1<float> Work(n);
    sgetri_(&n, 
            pInverse->data(), 
            &n, 
            pIPiv->data(),
            Work.data(),
            &n,
            &info);
    if (info!=0)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TSolverLapackLU::calculate_inverse: sgetri error %d\n",info);
      throw std::runtime_error(errormsg);
    }
    return pInverse;
  }


}

