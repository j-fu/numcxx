#ifndef TPRECON_JACOBI_HXX
#define TPRECON_JACOBI_HXX

#include "tsparsematrix.hxx"

namespace numcxx
{

    /// Bridge class for using umfpack as solver for vmatrix
    ///    
    /// UMFPACK is a GPL licensed sparse  matrix package by Tim Davis of Texas
    /// A&M university,  see [wikipedia](http://en.wikipedia.org/wiki/UMFPACK)
    /// and its [homepage](http://faculty.cse.tamu.edu/davis/suitesparse.html)
    template<typename T> 
    class TPreconJacobi: public  TLinSolver<T>
    {
        /// The corresponding matrix
        const std::shared_ptr< TSparseMatrix<T> > pMatrix;
        

    public:
        std::shared_ptr< TArray1<T> > pInvDiag;

        /// Create LU factorization class
        TPreconJacobi(const std::shared_ptr<TSparseMatrix<T>> pA);

        ~TPreconJacobi(){};
        /// Create LU factorization class

        static std::shared_ptr<TPreconJacobi<T>> create(const std::shared_ptr<TSparseMatrix<T>> pA);
        /// Perform actual computation of LU factorization
        void update();

        /// Solve LU factorized system
        void solve( TArray<T> & Sol,  const TArray<T> & Rhs) const;
    };
}

#include "tprecon-jacobi.ixx"

#endif

