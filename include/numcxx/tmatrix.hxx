#ifndef NUMCXX_TMATRIX_H
#define NUMCXX_TMATRIX_H

#include <cmath>
#include <cstdlib>

#include "tarray1.hxx"
#include "tarray2.hxx"

namespace numcxx
{
    /// Dense matrix class
    template<typename T> 
    class TMatrix: public TArray<T>, public TLinOperator<T>, public MatrixExpressionBase
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;
        using TArray<T>::operator[];
        
        /// Construct zero size array.
        TMatrix(): TArray<T>(){};
        
        /// Construct an empty matrix
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        TMatrix(index n): TArray<T>(n,n){};
        

        /// Construct matrix from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TMatrix(index n, T*data,std::function<void(T*p)> deleter):  TArray<T>(n,n,data,deleter){};
        
        /// Construct matrix from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method.
        /// \see TArray<T>#_datamanager
        TMatrix(index n, T*data, std::shared_ptr<void> datamanager): TArray<T>(n,n,data,datamanager){};
        

        /// Construct 2D Array from std::initializer list.
        TMatrix(const  std::initializer_list<std::initializer_list<T>> &il ): TArray<T>(il){};

        /// Copy constructor
        TMatrix(const TMatrix<T>& A):TArray<T>(A.shape(0),A.shape(1)){assign(*this,A);}


        
        /// Construct empty square matrix
        ///
        /// Mainly for access from python
        static std::shared_ptr<TMatrix <T> > create(index n) {  return std::make_shared<TMatrix<T>> (n);  }

        /// Construct matrix from std::initializer list.
        static std::shared_ptr<TMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il) { return std::make_shared<TMatrix<T>> (il);  }
        
        /// Create a copy of the matrix
        ///
        ///  \return Matrix of the same size with contents initialized to this
        std::shared_ptr<TMatrix <T> > copy() const  {  return std::make_shared<TMatrix<T>>(*this);  }
        
        /// Create a clone of the matrix
        ///
        ///  \return Matrix of the same size with empty contents.
        std::shared_ptr<TMatrix <T> > clone() const  {  return create(shape(0)); } 
        

        /// Matrix entry access for use in expression templates
        const T& xentry(const index i, const index j) const {return _data[_idx(i,j)];}
        
        /// Assignment operator
        TMatrix<T>&  operator=(const TMatrix<T> &expr) {return static_cast<TMatrix<T>&>(assign(*this,expr));}
        
        /// Assignment operator
        TMatrix<T>&  operator=(const T &expr) {return static_cast<TMatrix<T>&>(assign(*this,expr));}
        
        /// Apply matrix to vector
        ///
        /// \param u input
        /// \param v output: v=A*u
        void apply(const TArray<T> &u, TArray<T> &v) const;
        
        /// Calculate inverse of matrix
        std::shared_ptr<TMatrix<T>> calculate_inverse();
        
    private:
        using TArray<T>::_data;
        using TArray<T>::_idx;
        using TArray<T>::_check_square;
        
    };
    
}

#include "tmatrix.ixx"
#endif
        
