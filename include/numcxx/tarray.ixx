/// \file tarray.ixx
/// 
/// Inline metho definitions for numcxx::TArray
///
#include <cstdlib>
#include <cmath>

namespace numcxx
{

  template<typename T> 
  inline TArray<T>::TArray(const std::initializer_list<T> &il ):TArray(il.size())
  {
    index i=0;
    for (auto x = il.begin() ; x != il.end(); x++,i++) _data[i]= *x;
  }


  template <typename T> 
  inline TArray<T>::TArray(const  std::initializer_list<std::initializer_list<T>> &il):
    TArray(il.size(), il.begin()->size())
  {
    index i=0;
        
    for (auto jl = il.begin() ; jl != il.end(); jl++,i++)
    {
      index j=0;
      for (auto x = jl->begin() ; x != jl->end(); x++,j++) 
        _data[_idx(i,j)]= *x;
    }
  }

  template <typename T, typename EXPR,
            typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
  inline TArray<T>& assign(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
  {
    A.resize( expr.size() );
    T *data=A.data();
    for(index i=0; i<expr.size(); i++ ) data[i] = expr[i];
    return A;
  }

  template <typename T, typename VAL,
            typename= typename std::enable_if<!std::is_class<VAL>::value, VAL>::type>
  inline TArray<T>& assign(TArray<T>& A, const  VAL& a) 
  {
    T *data=A.data();
    for(index i=0; i<A.size(); i++ ) data[i] = a;
    return A;
  }

  template <typename T, typename EXPR,
            typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
  inline void xadd(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
  {
    A.resize( expr.size() );
    T *data=A.data();
    for(index i=0; i<expr.size(); i++ ) data[i] += expr[i];
  }


  template <typename T, typename VAL,
            typename= typename std::enable_if<!std::is_class<VAL>::value, VAL>::type>
  inline void xadd(TArray<T>& A, const  VAL& a) 
  {
    T *data=A.data();
    for(index i=0; i<A.size(); i++ )data[i] += a;
  }


  template <typename T, typename EXPR,
            typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
  inline void xsub(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
  {
    A.resize( expr.size() );
    T *data=A.data();
    for(index i=0; i<expr.size(); i++ ) data[i] -= expr[i];
  }


  template <typename T, typename VAL,
            typename= typename std::enable_if<!std::is_class<VAL>::value, VAL>::type>
  inline void xsub(TArray<T>& A, const  VAL& a) 
  {
    T *data=A.data();
    for(index i=0; i<A.size(); i++ )data[i] -= a;
  }


  template <typename T, typename EXPR,
            typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
  inline void xmul(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
  {
    A.resize( expr.size() );
    T *data=A.data();
    for(index i=0; i<expr.size(); i++ ) data[i] *= expr[i];
  }


  template <typename T, typename VAL,
            typename= typename std::enable_if<!std::is_class<VAL>::value, VAL>::type>
  inline void xmul(TArray<T>& A, const  VAL& a) 
  {
    T *data=A.data();
    for(index i=0; i<A.size(); i++ )data[i] *= a;
  }


  template <typename T, typename EXPR,
            typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
  inline void xdiv(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
  {
    A.resize( expr.size() );
    T *data=A.data();
    for(index i=0; i<expr.size(); i++ ) data[i] /= expr[i];
  }


  template <typename T, typename VAL,
            typename= typename std::enable_if<!std::is_class<VAL>::value, VAL>::type>
  inline void xdiv(TArray<T>& A, const  VAL& a) 
  {
    T *data=A.data();
    for(index i=0; i<A.size(); i++ )data[i] /= a;
  }

  template <typename T> 
  inline void TArray<T>::_check_bounds(index acc_dim, index acc_ndim, index acc_idx) const
  {
    if (acc_ndim!=_ndim) 
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TArray::_check_bounds: attempt of %uD access of %uD array",acc_ndim,_ndim);
      throw std::out_of_range(errormsg);
    }
    if ((acc_idx<0) || (acc_idx>=_shape[acc_dim]))
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TArray::_check_bounds: _shape[%u]=%u but i%u=%u",acc_dim,_shape[acc_dim],acc_dim,acc_idx);
      throw std::out_of_range(errormsg);
    }
  }  
    
  template <typename T> 
  inline void TArray<T>::_assert_square() const
  {
    if (_ndim!=2 || _shape[0]!=_shape[1])
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TArray::_assert_square: unexpected non-equal array dimensions\n");
      throw std::length_error(errormsg);
    }
  }

  template <typename T> 
  inline index TArray<T>::_idx(index i0) const
  {
#ifdef NUMCXX_CHECK_BOUNDS
    _check_bounds(0,1,i0);
#endif
    return i0;
  }

  template <typename T> 
  inline index TArray<T>::_idx(index i0,index i1)  const
  { 
#ifdef NUMCXX_CHECK_BOUNDS
    _check_bounds(0,2,i0);
    _check_bounds(1,2,i1);
#endif
    return i0*_shape[1]+i1;
  }
    
  template <typename T> 
  inline index TArray<T>::_idx(index i0,index i1, index i2)  const
  { 
#ifdef NUMCXX_CHECK_BOUNDS
    _check_bounds(0,3,i0);
    _check_bounds(1,3,i1);
    _check_bounds(2,3,i2);
#endif
    return (i0*_shape[0]+i1)*_shape[1]+i2;
  }
    
  template <typename T> 
  inline TArray<T>::TArray(index n0, index n1):
    _ndim(2),
    _shape{n0,n1},
    _size((size_t)n0*(size_t)n1),
    _data((T*)malloc(sizeof(T)*_size)),
    _deleter([](T*p){free(p);})
    {if (_data==nullptr)  throw std::runtime_error("numcxx: TArray::TArray(): Memory allocation failed"); };
    
  template <typename T> 
  inline TArray<T>::TArray(index n0):
    _ndim(1),
    _shape{n0},
    _size(n0),
    _data((T*)malloc(sizeof(T)*_size)),
    _deleter([](T*p){free(p);})
    {if (_data==nullptr)  throw std::runtime_error("numcxx: TArray::TArray(): Memory allocation failed"); };
    
  template <typename T> 
  inline TArray<T>::TArray(index n0, T*data,std::function<void(T*p)> deleter):
    _data(data),
    _ndim(1),
    _shape{n0},
    _size(n0),
    _deleter(deleter)
    {};
    
  template <typename T> 
  inline TArray<T>::TArray(index n0, T*data,std::shared_ptr<void> datamanager):
    TArray(n0,data,[](T*p){;}){_datamanager=datamanager;};
    
  template <typename T> 
  inline TArray<T>::TArray(index n0, index n1,T*data, std::function<void(T*p)> deleter):
    _data(data),
    _ndim(2),
    _shape{n0,n1},
    _size(n0*n1),
    _deleter(deleter)
    {};



    
  template <typename T> 
  inline TArray<T>::TArray(index n0, index n1, T*data,std::shared_ptr<void> datamanager):
    TArray(n0,n1,data, [](T*p){;})
  {_datamanager=datamanager;};
    
  template <typename T> 
  inline TArray<T>::~TArray()
  {
    if (_datamanager==nullptr)
    {
      _deleter(_data);
    }
    // otherwise we assume that datamaager takes care
    // of the data pointer
  };
    
  template <typename T> 
  inline TArray<T>::TArray():
    _data(0),
    _size(0),
    _ndim(1),
    _shape{0,0},
    _deleter([](T*p){;})
    {};   

  template <typename T> 
  inline void TArray<T>::_nullify()
  {
    _shape[0]=0;     
    _shape[1]=0;     
    _shape[2]=0;     
    _size=0;       
    _deleter=[](T*p){;};     
    _datamanager=nullptr; 
    _data=nullptr;
  }

  template <typename T> 
  inline void TArray<T>::_setshape(index shape0)
  {
    if (_ndim!=1)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TArray::resize: unable to set 1D shape for 2D array.\n");
      throw std::runtime_error(errormsg);
    }
    _shape[0]=shape0;   
    _shape[1]=0;
    _shape[2]=0;
    _size=shape0;       
  }
  

  template <typename T> 
  inline void TArray<T>::resize(size_t n)
  {
    if (_size==n) return;

    if (_ndim>1)
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TArray::resize: unable to resize 2D Array to 1D.\n");
      throw std::runtime_error(errormsg);
    }
        
    if (_datamanager==nullptr)
    {
      _deleter(_data);
      _data=(T*)malloc(sizeof(T)*n);
      if (_data==nullptr)  throw std::runtime_error("numcxx: TArray::resize(): Memory allocation failed"); 
      _deleter=[](T*p){free(p);};
      _size=n;
      _shape[0]=n;
      _shape[1]=0;
    }
    else
    {
      char errormsg[80];
      snprintf(errormsg,80,"numcxx::TArray::resize: unable to resize - data managed by different object.\n");
      throw std::runtime_error(errormsg);
    }
  }

    

  template <typename T> 
  inline void TArray<T>::operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B)  
  { 
    for(index i=0;i<A._size;i++) f(A._data[i],B._data[i]);
  }

  template <typename T> 
  inline void TArray<T>::operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C)  
  {
    for(index i=0;i<A._size;i++) f(A._data[i],B._data[i],C._data[i]);
  }

    
    
  template<typename T> 
  inline std::ostream & operator << (std::ostream & s, TArray<T> &A)
  {
    if (A.ndim()==1)
      for (index i=0;i<A.size();i++) s <<"[" << i << "]: " <<A(i) << std::endl << std::flush;
    else
    {
      s << "    ";
      for (index j=0;j<A.shape(1);j++) 
        s << "[" << j << "]     ";
      s<< "\n";
      for (index i=0;i<A.shape(0);i++) 
      {
        s << "[" << i << "]: ";
        for (index j=0;j<A.shape(1);j++) 
          s << A(i,j) << "   ";
        s<< "\n";
      }
    }
    return s;
  }


  template <typename T> 
  inline void TArray<T>::savetxt(std::ostream &s) const
  {
    if (ndim()==1)
      for (index i=0;i<size();i++) s << _data[i] << std::endl << std::flush;
    else
    {
      for (index i=0;i<shape(0);i++) 
      {
        for (index j=0;j<shape(1);j++) 
          s << _data[_idx(i,j)] << " ";
        s<< std::endl;
      }
      s << std::flush;
    }
  }

}
