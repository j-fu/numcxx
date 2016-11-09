namespace numcxx
{
    template<typename T> 
    inline TSolverLapackLU<T>::TSolverLapackLU(const std::shared_ptr<TMatrix<T>> a):
        TLinSolver<T>(),
        a(a),
        lu(a->clone()),
        ipiv(TArray1<int>::create(a->shape(0)))
    { update();;}
    
    template<typename T> 
    inline std::shared_ptr<TSolverLapackLU<T>> TSolverLapackLU<T>::create(const std::shared_ptr<TMatrix<T>> a)
    {
        return std::make_shared<TSolverLapackLU<T>>(a);
    }

    extern "C" 
    {
        void sgetrf_(int *n, int *m, float *a, int *lda, int* ipiv, int *info);
        void sgetrs_(char *trans,int *n, const int *nrhs, float*a, int* lda, int *ipiv , float *b, int *ldb, int *info );
        void dgetrf_(int *n, int *m, double *a, int *lda, int* ipiv, int *info);
        void dgetrs_(char *trans,int *n, const int *nrhs, double*a, int* lda, int *ipiv , double *b, int *ldb, int *info );
     }

    template<> 
    inline void TSolverLapackLU<double>::update()
    {
        int n=lu->shape(0);
        int info;
        *lu=*a;
        dgetrf_(&n,&n,lu->data(),&n,ipiv->data(),&info);
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
        int n=lu->shape(0);
        int one=1;
        int info;
        dgetrs_(trans,&n,&one,lu->data(),&n,ipiv->data(),sol.data(),&n,&info);
        if (info!=0)
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: dgetrs error %d\n",info);
            throw std::runtime_error(errormsg);
        }
    }


    template<> 
    inline void TSolverLapackLU<float>::update()
    {
        int n=lu->shape(0);
        int info;
        *lu=*a;
        sgetrf_(&n,&n,lu->data(),&n,ipiv->data(),&info);
        if (info!=0)
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: sgetrf error %d\n",info);
            throw std::runtime_error(errormsg);
        }

    }
    
    template<> 
    inline void TSolverLapackLU<float>::solve( TArray<float> & sol,  const TArray<float> & rhs) const
    {
        assign(sol,rhs);
        char trans[2]={'T','\0'};
        int n=lu->shape(0);
        int one=1;
        int info;
        sgetrs_(trans,&n,&one,lu->data(),&n,ipiv->data(),sol.data(),&n,&info);
        if (info!=0)
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TSolverLapackLU::update: sgetrs error %d\n",info);
            throw std::runtime_error(errormsg);
        }

    }

}

