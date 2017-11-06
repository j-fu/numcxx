#include "tsolver-lapacklu.hxx"
namespace numcxx
{


  extern "C" 
  {
    void dgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  double *ALPHA, const double *A, int *LDA,  double *B, int *LDB,  double *BETA,double * C, int *LDC);
    void sgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  float  *ALPHA, const float  *A, int *LDA,  float  *B, int *LDB,  float  *BETA,float  * C, int *LDC);
  }

  template <> 
  inline   void TMatrix<double>::apply(const TArray<double> &u, TArray<double> &v) const
  {
    v.resize(u.size());
    char transmat[2]={'T','\0'};
    char transvec[2]={'N','\0'};
    int n=shape(0);
    int ione=1;
    double done=1.0;
    double dzero=0.0;
    dgemm_(transmat,
           transvec,
           &n,&ione,&n,
           &done,
           _data,&n,
           u.data(),&n,
           &dzero,
           v.data(),&n);
  }


  template <> 
  inline   void TMatrix<float>::apply(const TArray<float> &u, TArray<float> &v) const
  {
    v.resize(u.size());
    char transmat[2]={'T','\0'};
    char transvec[2]={'N','\0'};
    int n=shape(0);
    int ione=1;
    float done=1.0;
    float dzero=0.0;
    sgemm_(transmat,
           transvec,
           &n,&ione,&n,
           &done,
           _data,&n,
           u.data(),&n,
           &dzero,
           v.data(),&n);
  }


  template <typename T>
  inline   void TMatrix<T>::apply(const TArray<T> &u, TArray<T> &v) const
  {
    for (index i=0;i<shape(1);i++)
    {   v[i]=0.0;
      for (index j=0;j<shape(0);j++)
        v[i]+=_data[_idx(i,j)]*u[j];
    }
  }

  template <typename T>
  inline std::shared_ptr<TMatrix<T>> TMatrix<T>::calculate_inverse()
  {
    auto pA=copy();  // essentially superfluous, but we seldomly use this
    auto pLU=TSolverLapackLU<T>::create(pA);
    return pLU->calculate_inverse();
  }
  
    
}

