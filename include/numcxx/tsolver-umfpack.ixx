#include <suitesparse/umfpack.h>
namespace numcxx
{
    
    /// Create LU factorization class
    template<typename T> 
    inline TSolverUMFPACK<T>::TSolverUMFPACK(const std::shared_ptr<TSparseMatrix<T>> pMatrix):
        pMatrix(pMatrix),
        Symbolic(nullptr),
        Numeric(nullptr)
    {}
    template<typename T> 
    inline TSolverUMFPACK<T>::~TSolverUMFPACK()
    {
        if (Symbolic!=nullptr) umfpack_di_free_symbolic(&Symbolic);
        if (Numeric!=nullptr) umfpack_di_free_numeric(&Numeric);
    }
    
    /// Create LU factorization class
    template<typename T> 
    inline std::shared_ptr<TSolverUMFPACK<T>> TSolverUMFPACK<T>::create(const std::shared_ptr<TSparseMatrix<T>> pA)
    {
        return std::make_shared<TSolverUMFPACK<T>>;
    }
    
    /// Perform actual computation of LU factorization
    template<typename T> 
    inline void TSolverUMFPACK<T>::update()
    {
        pMatrix->flush();
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
            if (Symbolic) umfpack_di_free_symbolic(&Symbolic),Symbolic=0;
            if (Numeric) umfpack_di_free_numeric(&Numeric),Numeric=0;
            status=umfpack_di_symbolic (n, n, pMatrix->pIA->data(), pMatrix->pJA->data(), pMatrix->pA->data(), &Symbolic, 0, 0);
            //      if (status>1) printf("Umfpack error %d\n",status);
        }
        if (Numeric!=nullptr) 
        {
            umfpack_di_free_numeric(&Numeric);
            Numeric=0;
        }
        status =umfpack_di_numeric (pMatrix->pIA->data(), pMatrix->pJA->data(), pMatrix->pA->data(), Symbolic, &Numeric, control, 0) ;
        //  if (status>1) printf("Umfpack error %d\n",status);
    }
    
    /// Solve LU factorized system
    template<typename T> 
    inline void TSolverUMFPACK<T>::solve( TArray<T> & Sol,  const TArray<T> & Rhs)
    {
        double *control=nullptr;
        int status;
        status=umfpack_di_solve (UMFPACK_At,pMatrix->pIA->data(), pMatrix->pJA->data(), pMatrix->pA->data(), Sol.data(), Rhs.data(), Numeric, control, 0 ) ;
//        if (status>1) printf("Umfpack error %d\n",status);
        
    }


}
