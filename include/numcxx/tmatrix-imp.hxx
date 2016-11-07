namespace numcxx
{
    template <typename T> 
    inline TMatrix<T>::TMatrix(): TArray2<T>(){};
    
    template <typename T> 
    inline TMatrix<T>::TMatrix(index n): TArray2<T>(n,n){};
    
    template <typename T> 
    inline TMatrix<T>::TMatrix(index n, T*data,std::function<void(T*p)> deleter): 
        TArray2<T>(n,n,data,deleter){};
    
    template <typename T> 
    inline TMatrix<T>::TMatrix(index n, T*data, std::shared_ptr<void> datamanager): 
        TArray2<T>(n,n,data,datamanager){};

    template <typename T> 
    inline TMatrix<T>::TMatrix(const  std::initializer_list<std::initializer_list<T>> &il ): TArray2<T>(il){};

    template <typename T> 
    inline std::shared_ptr<TMatrix <T> > TMatrix<T>::create(index n)
    {
        return std::make_shared<TMatrix<T>> (n);
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
        *x=*this;
        return x;
    }

    template <typename T> 
    inline   std::shared_ptr<TMatrix <T> > TMatrix<T>::clone() const
    { 
        return create(shape(0));
    } 


    extern "C" 
    {
        void dgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  double *ALPHA, const double *A, int *LDA,  double *B, int *LDB,  double *BETA,double * C, int *LDC);
        void sgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  float  *ALPHA, const float  *A, int *LDA,  float  *B, int *LDB,  float  *BETA,float  * C, int *LDC);
    }

    template <> 
    inline   void TMatrix<double>::apply(const TArray1<double> &u, TArray1<double> &v) const
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
    inline   void TMatrix<float>::apply(const TArray1<float> &u, TArray1<float> &v) const
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


    template <typename T>
    inline   void TMatrix<T>::apply(const TArray1<T> &u, TArray1<T> &v) const
    {
        for (index i=0;i<shape(1);i++)
        {   v[i]=0.0;
            for (index j=0;j<shape(0);j++)
                v[i]+=_data[_idx(i,j)]*u[j];
        }
    }
    
}

