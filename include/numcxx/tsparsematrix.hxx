#ifndef NUMCXX_TSPARSEMATRIX_H
#define NUMCXX_TSPARSEMATRIX_H

#include <vector>
#include "tarray1.hxx"
#include "tmatrix.hxx"

namespace numcxx
{
    
  template<typename T>   class TSolverUMFPACK;
  template<typename T>   class TPreconJacobi;

  /// Sparse matrix class using CRS storage scheme
  template<typename T> 
  class TSparseMatrix: public TLinOperator<T>, public SparseMatrixExpressionBase
  {
    friend class TSolverUMFPACK<T>;
    friend class TPreconJacobi<T>;
  public:
    typedef T value_type;

    TSparseMatrix():n(0) {};

    /// Create an empty sparse matrix representing
    /// an  n1 x n2 system of linear equations
    /// (currently, n1 must be equal to n2)
    TSparseMatrix(index n1, index n2);


    /// Static wrapper around corresponding constructor
    static std::shared_ptr<TSparseMatrix <T> > create(index n1, index n2) {return std::make_shared<TSparseMatrix<T>>(n1,n2);}

    /// Create an empty sparse matrix representing
    /// an  n x n system of linear equations
    TSparseMatrix(const  std::initializer_list<std::initializer_list<T>> &il);

    /// Static wrapper around corresponding constructor
    static std::shared_ptr<TSparseMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);


    /// Access operator. 
    ///
    /// This operator accesses the sparse matrix element at position
    /// i,j. If it is not existent, it is created. All freshly created
    /// elements are collected in an extension list which is joined with
    /// the main data structure using the flush() method.
    T& operator()(int i, int j);
        
    /// Re-create the internal data structure in order to
    /// accomodated all newly created elements.
    void  flush();


    void  clear() {  flush(); (*pIA)=0.0;}

    /// Apply sparse matrix to vector
    void apply(const TArray<T> &U,TArray<T> &V ) const;

    /// Create a dense matrix from the sparse matrix
    std::shared_ptr<TMatrix<T>> copy_as_dense();

    /// Return the shape od the matrix
    index shape(int idim) {return n;}
        
    /// Copy constructor is deleted
    TSparseMatrix(const TSparseMatrix<T>& A)=delete;


    // std::shared_ptr<TSparseMatrix <T> > copy() const;
    // std::shared_ptr<TSparseMatrix <T> > clone() const;
    // T xentry(const index i, const index j) const {return _data[_idx(i,j)];}
    // template <typename VAL>
    // TSparseMatrix<T>&  operator=(const VAL  &expr)  { assign(*this,expr); return *this;}
    // TSparseMatrix<T>&  operator=(const TSparseMatrix<T> &expr) { assign(*this,expr); return *this;}
    // void apply(const TArray1<T> &u, TArray1<T> &v) const;

    /// Row pointers
    std::shared_ptr<TArray1<int>> pIA;
        
    /// Column indices
    std::shared_ptr<TArray1<int>> pJA;  
        
    /// Entries
    std::shared_ptr<TArray1 <T> > pA; 

    /// Calculate inverse of sparse matrix whuch is dense
    std::shared_ptr<TMatrix<T>> calculate_inverse();

  private:



    /// Check if pattern has changed after last solver
    /// update
    bool pattern_changed(){return _pattern_changed;}

    void pattern_changed(bool chg) {_pattern_changed=chg;};

        
    bool _pattern_changed=false;
    const index n;
    int maxrow=0;

    bool _first_flush_done=false;
    bool empty() { return !_first_flush_done;}
    /// Row entry for flush method
    class RowEntry;

    /// Extension data
    class Extension;

    std::shared_ptr<Extension>  pExt=nullptr;
        
  };
    
}

#include "tsparsematrix.ixx"
#endif
