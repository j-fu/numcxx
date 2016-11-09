namespace numcxx
{
    
    /// Create LU factorization class
    template<typename T> 
    inline TPreconJacobi<T>::TPreconJacobi(const std::shared_ptr<TSparseMatrix<T>> pMatrix):
        pMatrix(pMatrix)
    {
        pInvDiag=std::make_shared<TArray1<T>>(pMatrix->shape(0));
        update();
    }

    
    /// Create LU factorization class
    template<typename T> 
    inline std::shared_ptr<TPreconJacobi<T>> TPreconJacobi<T>::create(const std::shared_ptr<TSparseMatrix<T>> pA)
    {
        return std::make_shared<TPreconJacobi<T>>(pA);
    }
    
    /// Perform actual computation of LU factorization
    template<typename T> 
    inline void TPreconJacobi<T>::update()
    {
        auto &M=*pMatrix;
        int n=M.shape(0);
        auto &InvDiag=*pInvDiag;
        for (int i=0;i<n;i++)
            InvDiag(i)=1.0/M(i,i);

        pMatrix->pattern_changed(false);
    }
    
    /// Solve LU factorized system
    template<typename T> 
    inline void TPreconJacobi<T>::solve( TArray<T> & Sol,  const TArray<T> & Rhs) const
    {
        Sol.resize(Rhs.size());
        int n=pMatrix->shape(0);
        auto &InvDiag=*pInvDiag;
        for (int i=0;i<n;i++)
            Sol(i)=Rhs(i)*InvDiag(i);
    }


}
