/// \file tarray.hxx
/// 
/// Header for numcxx::TArray
///
#ifndef NUMCXX_TARRAY_H
#define NUMCXX_TARRAY_H
#include <ostream>
#include <typeinfo>
#include <memory>
#include <stdexcept> 
#include "expression.ixx"


namespace numcxx
{
  template<typename T> class TArray;



  template<typename T>
  inline std::ostream & operator << (std::ostream & s, TArray<T> &A);

  /// TArray is the common template base class for arrays and dense matrices
  /// of the numcxx project.
  template<typename T> class TArray
  {
  public:

    typedef T value_type;

    /// Access operator for 1D arrays.
    ///
    /// \param i0  Index of element to be accessed.
    /// \return    Reference to element to be accessed.
    T & operator()(const index i0)  { return _data[_idx(i0)];};
    const T & operator()(const index i0) const  { return _data[_idx(i0)];};

    /// Access operator for 2D arrays.
    ///
    /// \param i0  Row index of element to be accessed.
    /// \param i0  Column index of element to be accessed.
    /// \return    Reference to element to be accessed.
    T & operator()(const index i0, const index i1)  { return _data[_idx(i0,i1)];};
    const T & operator()(const index i0, const index i1) const { return _data[_idx(i0,i1)];};
    

    /// Obtain tensor dimension of array.
    /// 
    /// Tensor dimension is 1 for vectors, 2 for matrices.
    /// \return Dimension.
    index ndim() const {return _ndim;}

    /// Obtain size of array.
    /// 
    /// This ist the overall number of elements in the array
    /// \return Size.
    size_t size() const { return _size;}

    /// Obtain shape of array for given dimension.
    ///
    /// For 1D arrays, ``shape(0)``  is equivalent to the size
    /// For 2D arrays, ``shape(0)`` is the number of rows
    /// and ``shape(1)`` the number of columns. This corresponds
    /// to the "row major" storage format.
    /// \param dim Tensor dimension.
    /// \return Number of elements in given dimension.
    index shape(const index dim)  const {return _shape[dim];}

    /// Add value to all elements.
    ///
    /// \param a Summand for each element.
    template <typename VAL>
    void operator+=(const VAL & a) {xadd(*this,a);}

    /// Subtract value from all elements.
    ///
    /// \param a Value to be subracted from each element.
    template <typename VAL>
    void operator-=(const VAL& a) {xsub(*this,a);}


    /// Multiply all elements by value
    ///
    /// \param a Multiplicator for each element.
    template <typename VAL>
    void operator*=(const VAL& a) {xmul(*this,a);}


    /// Divide each element by value
    ///
    /// \param a Divisor for each element.
    template <typename VAL>
    void operator/=(const VAL & a) {xdiv(*this,a);}



    /// Binary operation on arrays.
    ///
    /// \param f Function performing operation for each index.
    /// \param A First array argument.
    /// \param B Second array argument.
    static void operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B);

    /// Ternary operation on arrays.
    ///
    /// \param f Function performing operation for each index.
    /// \param A First array argument.
    /// \param B Second array argument.
    /// \param C Third array argument.
    static void operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C);


    /// Alternative access operator for 1D arrays
    T & operator[](const index i) { return _data[i];}

    /// Const reference to entry for use in expression templates
    const T & operator[](const index i) const  { return _data[i];};
        
    /// Expression template compatible assignment operator
    template <typename VAL>
    TArray<T>&  operator=(const VAL  &expr)  {return assign(*this,expr);}

    /// Obtain C-pointer of data array.
    /// 
    /// \return Address of C-Array managed by the class which holds
    ///         the data
    T*data() const { return _data;}


    /// Resize array
    void resize(size_t n);


    /// Copy constructor is deleted
    TArray(const TArray<T>& A)=delete;


  private:

    /// Tensor dimension.
    const index _ndim;

    /// Size of array.
    size_t _size;
        
    /// Shape vector
    index _shape[3]={0,0,0};


    /// Bounds checker.
    /// 
    /// Bounds check is enabled if the code is compiled with
    /// ``-DNUMCXX_CHECK_BOUNDS``.
    ///
    /// On error it throws ``std::out_of_range`` exception.
    /// \param acc_dim  Dimension to be checked.
    /// \param acc_ndim Tensor dimension to be checked.
    /// \param acc_dim  Index in dimension \p acc_dim to be checked.
    void _check_bounds(index acc_dim, index acc_ndim, index acc_idx) const;


  protected:

    /// Deleter method.
    /// 
    /// This is the  proper method to be used to  destroy the data
    /// pointer in the array if data manager is null. Depending on
    /// the way  it was constructed,  it may do  nothing, ``free``
    /// the memory, ``delete[]`` the memory, or something else.
    std::function<void(T*p)> _deleter=nullptr;

    /// Data manager.
    /// 
    /// Smart pointer to some other object managing the data pointer.
    /// If it is not nullptr, the deleter is not called and the 
    /// memory corresponding to the data pointer is freed when the
    /// object behind the data manager is destroyed.
    ///
    /// An example in case is the use of ``shared_ptr<vector> v``
    /// as datamanager and  ``v->data()`` as data pointer
    std::shared_ptr<void>_datamanager =nullptr;

    /// Data pointer.
    T* _data=nullptr;

       
    /// Check if all shapes are the same.
    ///
    /// Throw an exception on error
    void _assert_square() const;
        
    /// 1D Array index calculation with optional bounds check.
    index _idx(index i0) const;

    /// 2D Array index calculation with optional bounds check.
    index _idx(index i0,index i1)  const;

    /// 3D Array index calculation with optional bounds check.
    index _idx(index i0,index i1,index i2)  const;

    /// Construct an zero length 1D array.
    TArray();

    ///  Nullify contents of array (for move constructors)
    void _nullify();

    ///  Set shape of 1D array (for move constructors)
    void _setshape(index shape0);

    /// Construct an empty 1D array of length n0
    /// \param n0 Size.
    TArray(index n0);

    /// Construct 1D Array from std::initializer list.
    TArray(const std::initializer_list<T> &il );


    /// Construct an 1D array from data pointer
    ///
    /// \param n0 Size.
    /// \param data Pointer to data.
    /// \param deleter Deleter method, \see TArray<T>#_deleter
    TArray(index n0, T*data, std::function<void(T*p)> deleter);

    /// Construct an 1D array from data pointer
    ///
    /// \param n0 Size.
    /// \param data Pointer to data.
    /// \param deleter Deleter method.
    /// \see TArray<T>#_datamanager
    TArray(index n0, T*data, std::shared_ptr<void> datamanager);

    /// Construct an empty 2D array.
    ///
    /// \param n0 Number of rows
    /// \param n1 Number of columns
    TArray(index n0, index n1);

    /// Construct 2D Array from std::initializer list.
    TArray(const  std::initializer_list<std::initializer_list<T>> &il );


    /// Construct a 2D array from data pointer
    ///
    /// \param n0 Number of rows
    /// \param n1 Number of columns
    /// \param data Pointer to data.
    /// \param deleter Deleter method, \see TArray<T>#_deleter
    TArray(index n0, index n1, T*data,std::function<void(T*p)> deleter);

    /// Construct a 2D array from data pointer
    ///
    /// \param n0 Number of rows
    /// \param n1 Number of columns
    /// \param data Pointer to data.
    /// \param deleter Deleter method.
    /// \see TArray<T>#_datamanager
    TArray(index n0, index n1, T*data,std::shared_ptr<void> datamanager);
        

    /// Destructor.
    ~TArray();

    /// Print contents of array.
    friend std::ostream & operator<< <T>(std::ostream & s, TArray<T> &A);


  };



  /// Base class for linear solvers and preconditioners
  template<typename T> class TLinSolver
  {
  public:
    TLinSolver(){};
    virtual void solve( TArray<T> & sol,  const TArray<T> & rhs) const {};
    virtual void update(void){};
  };

  /// Base class for linear operators (matrices, sparse matrices)
  template<typename T> class TLinOperator
  {
  public:
    TLinOperator(){};
    virtual void apply( const TArray<T> & sol,   TArray<T> & rhs) const {};
  };
    
    
}


#include "tarray.ixx"
#endif
