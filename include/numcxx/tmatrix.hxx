#ifndef NUMCXX_TMATRIX_H
#define NUMCXX_TMATRIX_H

#include <cmath>
#include <cstdlib>

#include "tarray.hxx"

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
    TMatrix(): TArray<T>(){_assert_square();};
        
    /// Construct an empty matrix
    ///
    /// \param n0 Number of rows
    /// \param n1 Number of columns
    TMatrix(index n1, index n2): TArray<T>(n1,n2){_assert_square();};
        

    /// Construct matrix from data pointer
    ///
    /// \param n0 Number of rows
    /// \param n1 Number of columns
    /// \param data Pointer to data.
    /// \param deleter Deleter method, \see TArray<T>#_deleter
    TMatrix(index n1, index n2, T*data,std::function<void(T*p)> deleter):  TArray<T>(n1,n2,data,deleter){_assert_square();};
        
    /// Construct matrix from data pointer
    ///
    /// \param n0 Number of rows
    /// \param n1 Number of columns
    /// \param data Pointer to data.
    /// \param deleter Deleter method.
    /// \see TArray<T>#_datamanager
    TMatrix(index n1, index n2, T*data, std::shared_ptr<void> datamanager): TArray<T>(n1,n2,data,datamanager){_assert_square();};
        

    /// Construct 2D Array from std::initializer list.
    TMatrix(const  std::initializer_list<std::initializer_list<T>> &il ): TArray<T>(il){_assert_square();};

    /// Copy constructor
    TMatrix(const TMatrix<T>& A):TArray<T>(A.shape(0),A.shape(1)){assign(*this,A);}

        
    /// Construct empty square matrix
    ///
    /// Mainly for access from python
    static std::shared_ptr<TMatrix <T> > create(index n1, index n2) {  return std::make_shared<TMatrix<T>> (n1,n2);  }

    /// Construct matrix from std::initializer list.
    static std::shared_ptr<TMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il) { return std::make_shared<TMatrix<T>> (il);  }
        
    /// Create a copy of the matrix
    ///
    ///  \return Matrix of the same size with contents initialized to this
    std::shared_ptr<TMatrix <T> > copy() const  {  return std::make_shared<TMatrix<T>>(*this);  }
        
    /// Create a clone of the matrix
    ///
    ///  \return Matrix of the same size with empty contents.
    std::shared_ptr<TMatrix <T> > clone() const  {  return create(shape(0),shape(1)); } 
        
    /// Element read access.
    /// 
    /// \param i0  Row index of element to be accessed.
    /// \param i1  Column index of element to be accessed.
    /// \return Value of element at (i0,i1)
    T item(index i0,index i1)  { return _data[_idx(i0,i1)];};


    /// Element write access.
    /// 
    /// \param i0  Row index of element to be accessed.
    /// \param i1  Column index of element to be accessed.
    /// \param x value to be copied to element at index.
    void itemset(index i0, index i1, T x)  { _data[_idx(i0,i1)]=x;};


    /// Matrix entry access for use in expression templates
    const T& xentry(const index i, const index j) const {return _data[_idx(i,j)];}

    /// Getter routine for access from python.
    /// 
    /// This access is rather expensive, as it constructs
    /// a smart pointer to the row. 
    /// \param i0  row index of element to be accessed.
    /// \return Smart pointer to i0-th row. 
    std::shared_ptr<TArray1 <T> > const __getitem__(index i0){ return std::shared_ptr<TArray1<T>>(new TArray1<T>(shape(1), &_data[_idx(i0,0)], [](T*p){;}));}


        
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

    bool is_matrix(){return true;}
        
  private:
    using TArray<T>::_data;
    using TArray<T>::_idx;
    using TArray<T>::_assert_square;
        
  };
    
}

#include "tmatrix.ixx"
#endif
        
