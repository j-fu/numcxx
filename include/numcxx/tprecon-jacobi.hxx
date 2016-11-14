#ifndef TPRECON_JACOBI_HXX
#define TPRECON_JACOBI_HXX

#include "tsparsematrix.hxx"

namespace numcxx
{

    /// Jacobi preconditioner class
    template<typename T> 
    class TPreconJacobi: public  TLinSolver<T>
    {
        /// The corresponding matrix
        const std::shared_ptr< TSparseMatrix<T> > pMatrix;
        

    public:
        std::shared_ptr< TArray1<T> > pInvDiag;

        /// Create Preconditioner
        TPreconJacobi(const std::shared_ptr<TSparseMatrix<T>> pA);

        ~TPreconJacobi(){};

        /// Create preconditioner
        static std::shared_ptr<TPreconJacobi<T>> create(const std::shared_ptr<TSparseMatrix<T>> pA);

        /// Perform actual computation preconditioner
        void update();

        /// Solve preconditioning system
        void solve( TArray<T> & Sol,  const TArray<T> & Rhs) const;
    };
}

#include "tprecon-jacobi.ixx"

#endif

