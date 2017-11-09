/// \file tarray1.hxx
/// 
/// Header for numcxx::TArray1
///
#ifndef NUMCXX_TARRAY1_H
#define NUMCXX_TARRAY1_H
#include <vector>
#include "tarray.hxx"

namespace numcxx
{

  ///
  /// One dimensional array class
  /// 
  /// Instances of this class can be created in various ways. The preferred
  /// construction of empty array goes like this:
  /// ````
  /// numcxx::TArray1<double> A(n);
  /// std::shared_ptr<numcxx::TArray1<double>> pA=numcxx::TArray1<double>::create(n)
  /// ````
  /// 
  /// As a derived class from numcxx::TArray<T>  it is merely a facade to the content
  /// in the base class. 
  /// 
  /// An alias numcxx::DArray1 for numcxx::TArray1<double> is available from numcxx.h.
  ///
  /// For this class, expression templates are defined in expression.ixx which allow
  /// to use these arrays in standard linear algebra expressions.
  ///
  template<typename T> class TArray1: public TArray<T>, public ExpressionBase
  {
  public:
    using TArray<T>::operator[];
    using TArray<T>::size;
    using TArray<T>::resize;
    using TArray<T>::operator=;


        

    /// Construct an empty 1D array.
    ///
    /// After construction, the entry values are unitialized
    /// \param n0 Size.
    TArray1(index n):TArray<T>(n){};

    /// Construct 1D Array from std::initializer list.
    /// 
    /// This is the preferable way for small arrays. It allows to write e.g.
    /// ````
    /// TArray1<double> A={1,2,3};
    /// ````
    TArray1(const std::initializer_list<T> &il ):TArray<T>(il){};

    /// Construct an 1D array from  plain old C Array
    ///
    /// \param n0 Size.
    /// \param data Pointer to data.
    /// \param deleter Deleter method, \see TArray<T>#_deleter
    TArray1(index n0, T*data, std::function<void(T*p)> deleter):TArray<T>(n0,data,deleter){};

    /// Construct an 1D array from  smartpointer managed array
    ///
    /// \param n0 Size.
    /// \param data Pointer to data.
    /// \param datamanager Smartpointer to object managing data.
    /// \see TArray<T>#_datamanager
    TArray1(index n0, T*data, std::shared_ptr<void> datamanager):TArray<T>(n0,data,datamanager){}; 


    /// Construct 1D Array from smart pointer to std::vector.
    ///
    ///  \param v Smart pointer to vector. A copy (increasing refcount)
    /// is stored as datamanager in created object.
    TArray1(std::shared_ptr<std::vector<T>> v):TArray<T>(v->size(),v->data(),v){};

    ///
    /// Construct smart pointer empty 1D Array
    ///
    static std::shared_ptr<TArray1 <T> > create(index n1) { return std::make_shared<TArray1 <T> >(n1);}


    /// Construct 1D Array from std::initializer list.
    ///
    /// Static factory method to be used in place of ``std::make_shared``
    /// which has problems to detect the type of the ``initializer_list``.
    static std::shared_ptr<TArray1 <T> > create(const std::initializer_list<T> il){return std::make_shared<TArray1 <T> >(il);}

    /// Create a copy
    ///
    std::shared_ptr<TArray1 <T> > copy() const {return std::make_shared<TArray1 <T>>(*this); }


    /// Create a clone of the array.
    ///
    ///  \return Array of the same size with empty contents.
    std::shared_ptr<TArray1 <T> > clone() const     {  return create(size()); }



    /// Assignment operator
    TArray1<T>&  operator=(const TArray1<T> &expr) {return static_cast<TArray1<T>&>(assign(*this,expr));}

    /// Copy constructor
    TArray1(const TArray1<T>& A):TArray1<T>(A.shape(0)){assign(*this,A);}

    /// Move constructor
    TArray1(TArray1<T> && A):TArray1<T>(A.shape(0), A._data, [](T*p){;})
    {
      TArray<T>::_deleter=std::move(A._deleter);
      TArray<T>::_datamanager=A._datamanager;  
      A._nullify();
    }

    /// Move assignment
    TArray1<T>& operator=(TArray1<T> && A){
      if ( TArray<T>::_datamanager==nullptr)   TArray<T>::_deleter(_data);
      TArray<T>::_setshape(A.shape(0));
      TArray<T>::_data=A._data; 
      TArray<T>::_deleter=A._deleter;
      TArray<T>::_datamanager=A._datamanager;  
      A._nullify();
      return *this;
    }
        
    /// Copy constructor from expression
    template <typename EXPR, typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
    TArray1(const EXPR& A):TArray1<T>(A.size()){assign(*this,A);}

    /// Element read access.
    /// 
    /// \param i0  Index of element to be accessed.
    /// \return Value of element at i0
    T item(index i0) const  { return _data[_idx(i0)];};


    /// Element write access.
    /// 
    /// \param i0  Index of element to be accessed.
    /// \param x value to be copied to element at index.
    void itemset(index i0, T x) { _data[_idx(i0)]=x;};


    /// Element read access.
    /// 
    /// Getter routine for access from python.
    /// \param i0  Index of element to be accessed.
    /// \return Value of element at i0
    T __getitem__(index i0) const  { return _data[_idx(i0)]; };


    /// Element write access.
    /// 
    /// Setter routine for access from python.
    /// \param i0  Index of element to be accessed.
    /// \param x value to be copied to element at index.
    void __setitem__(index i0,T x)  { _data[_idx(i0)]=x; };

    
    bool is_matrix()  {return false;};

    /// Default constructor.
    TArray1():TArray<T>(){};

  private:
    using TArray<T>::_data;
    using TArray<T>::_idx;
       
  };
}
#endif
