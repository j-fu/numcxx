
namespace numcxx
{
    template <typename T> 
    inline TMatrix<T>::TMatrix(): TArray2<T>(){};
    
    template <typename T> 
    inline TMatrix<T>::TMatrix(index n0, index n1): TArray2<T>(n0,n1){};
    
    template <typename T> 
    inline TMatrix<T>::TMatrix(index n0, index n1, T*data,std::function<void(T*p)> deleter): 
        TArray2<T>(n0,n1,data,deleter){};
    
    template <typename T> 
    inline TMatrix<T>::TMatrix(index n0,index n1, T*data, std::shared_ptr<void> datamanager): 
        TArray2<T>(n0,n1,data,datamanager){};

    template <typename T> 
    inline TMatrix<T>::TMatrix(const  std::initializer_list<std::initializer_list<T>> &il ): TArray2<T>(il){};

    template <typename T> 
    inline std::shared_ptr<TMatrix <T> > TMatrix<T>::create(index n0,index n1)
    {
        return std::make_shared<TMatrix<T>> (n0,n1);
    }
    
    template <typename T> 
    inline  std::shared_ptr<TMatrix <T> > TMatrix<T>::create(const  std::initializer_list<std::initializer_list<T>> &il)
    {
        return std::make_shared<TMatrix<T>> (il);
    }
    
    template <typename T> 
    inline  std::shared_ptr<TMatrix <T> > TMatrix<T>::copy() const
    {
        auto x=create(shape(0),shape(1));
        x->fill(*this);
        return x;
    }

    template <typename T> 
    inline   std::shared_ptr<TMatrix <T> > TMatrix<T>::clone() const
    { 
        return create(shape(0),shape(1));
    } 

    
    
    

    template <typename T> 
    inline   std::tuple<std::shared_ptr<TMatrix<T>>, std::shared_ptr<TArray1<int>>> TMatrix<T>::lu_decomp() const
    {
        _check_square();
        auto lu=this->copy();
        auto ipiv=TArray1<int>::create(shape(0));
        lu_decomp(*lu,*ipiv);
        return std::make_tuple(lu,ipiv);
    }

    template <typename T> 
    inline   std::shared_ptr<TArray1<T>> TMatrix<T>::lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,const TArray1<T> &rhs)
    {
        auto sol=rhs.clone();
        lu_solve(lu,ipiv,*sol,rhs);
        return sol;
    }


    template <typename T> 
    inline   void TMatrix<T>::solve(TArray1<T> &sol, const TArray1<T> &rhs) const
    {
        std::shared_ptr< TMatrix<T> > lu;
        std::shared_ptr< TArray1<int> > ipiv;
        std::tie(lu,ipiv)=lu_decomp();
        lu_solve(*lu,*ipiv,sol,rhs);
    }
    
    template <typename T> 
    inline   std::shared_ptr<TArray1<T>> TMatrix<T>::solve(const TArray1<T> &rhs) const
    {
        auto sol=rhs.clone();
        solve(*sol,rhs);
        return sol;
    }
    
    
    
    template <typename T> 
    inline   std::shared_ptr<TArray1<T>> TMatrix<T>::apply(const TArray1<T> &u)
    {
        auto v=u.clone();
        apply(u,*v);
        return v;
    }
    
    
    
    template <typename T> 
    inline void TMatrix<T>::inverse(TMatrix<T>& inv) const
    {
        std::shared_ptr< TMatrix<T> > lu;
        std::shared_ptr< TArray1<int> > ipiv;
        std::tie(lu,ipiv)=lu_decomp();
        inverse(*lu,*ipiv,inv);
    }
    
    
    template <typename T> 
    inline std::shared_ptr<TMatrix<T>> TMatrix<T>::inverse(const TMatrix<T> &lu, const TArray1<int> &ipiv)
    {
        auto inv=lu.clone();
        inverse(lu,ipiv,*inv);
        return inv;
    }

    template <typename T> 
    inline std::shared_ptr<TMatrix<T>> TMatrix<T>::inverse() const
    {
        auto inv=clone();
        inverse(*inv);
        return inv;
    }


    /////////////////////////////////////////////////////////////////////////
    // LAPACK INTERFACE FOR DOUBLE

    extern "C"
    {
        extern void dgetrf_(int *n, int *m, double *a, int *lda, int* ipiv, int *info);
        extern void dgetrs_(char *trans,int *n, const int *nrhs, double*a, int* lda, int *ipiv , double *b, int *ldb, int *info );
        extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  double *ALPHA, const double *A, int *LDA,  double *B, int *LDB,  double *BETA,double * C, int *LDC);
        extern void dgetri_(int *n, double*a, int* lda, int *ipiv , double *work, int *lwork, int *info );
    }

    template <> 
    inline   void TMatrix<double>::lu_decomp(TMatrix<double> &lu, TArray1<int> &ipiv) 
    {
        int n=lu.shape(0);
        int info;
        dgetrf_(&n,&n,lu.data(),&n,ipiv.data(),&info);
    }

    template <> 
    inline   void TMatrix<double>::lu_solve(TMatrix<double> &lu, TArray1<int> &ipiv,TArray1<double> &sol, const TArray1<double> &rhs)
    {
        sol.fill(rhs);
        char trans[2]={'T','\0'};
        int n=lu.shape(0);
        int one=1;
        int info;
        dgetrs_(trans,&n,&one,lu.data(),&n,ipiv.data(),sol.data(),&n,&info);
    }

    template <> 
    inline   void TMatrix<double>::apply(const TArray1<double> &u, TArray1<double> &v)
    {
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
    inline void TMatrix<double>::inverse(TMatrix<double> &lu, TArray1<int> &ipiv, TMatrix<double>& inv)
    {
        int n=lu.shape(0);
        int info;
        TArray1<double> work(n);
        inv.fill(lu);
        dgetri_(&n, 
                inv.data(), 
                &n, 
                ipiv.data(),
                work.data(),
                &n,
                &info);
    }
    
    /////////////////////////////////////////////////////////////////////////
    // LAPACK INTERFACE FOR FLOAT
    extern "C"
    {
        extern void sgetrf_(int *n, int *m, float *a, int *lda, int* ipiv, int *info);
        extern void sgetrs_(char *trans,int *n, const int *nrhs, float*a, int* lda, int *ipiv , float *b, int *ldb, int *info );
        extern void sgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  float *ALPHA, const float *A, int *LDA,  float *B, int *LDB,  float *BETA,float * C, int *LDC);
        extern void sgetri_(int *n, float*a, int* lda, int *ipiv , float *work, int *lwork, int *info );
    }

    template <> 
    inline   void TMatrix<float>::lu_decomp(TMatrix<float> &lu, TArray1<int> &ipiv) 
    {
        int n=lu.shape(0);
        int info;
        sgetrf_(&n,&n,lu.data(),&n,ipiv.data(),&info);
    }

    template <> 
    inline   void TMatrix<float>::lu_solve(TMatrix<float> &lu, TArray1<int> &ipiv,TArray1<float> &sol, const TArray1<float> &rhs)
    {
        sol.fill(rhs);
        char trans[2]={'T','\0'};
        int n=lu.shape(0);
        int one=1;
        int info;
        sgetrs_(trans,&n,&one,lu.data(),&n,ipiv.data(),sol.data(),&n,&info);
    }

    template <> 
    inline   void TMatrix<float>::apply(const TArray1<float> &u, TArray1<float> &v)
    {
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

    template <> 
    inline void TMatrix<float>::inverse(TMatrix<float> &lu, TArray1<int> &ipiv, TMatrix<float>& inv)
    {
        int n=lu.shape(0);
        int info;
        TArray1<float> work(n);
        inv.fill(lu);
        sgetri_(&n, 
                inv.data(), 
                &n, 
                ipiv.data(),
                work.data(),
                &n,
                &info);
    }



}
