/////////////////////////////////////////////////
// Inline methods of TArray<T>
#include <cstdlib>
#include <cmath>

namespace numcxx
{


    template <typename T> 
    inline        T*TArray<T>::data() const { return _data;}
    template <typename T> 
    inline    index TArray<T>::ndim() const {return _ndim;}
    template <typename T> 
    inline    index TArray<T>::size() const {return _size;}
    template <typename T> 
    inline    index TArray<T>::shape(const index dim)  const {return _shape[dim];}


    template <typename T> 
    inline void TArray<T>::fill(T x)
    {
        for (index i=0;i<_size;i++) _data[i]=x;
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
    inline void TArray<T>::_check_square() const
    {
        if (_ndim!=2 || _shape[0]!=_shape[1])
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TArray::_check_square: unexpected non-equal array dimensions\n");
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
        _size(n0*n1),
        _data((T*)malloc(sizeof(T)*_size)),
        _deleter([](T*p){free(p);})
        {};
    
    template <typename T> 
        inline TArray<T>::TArray(index n0):
        _ndim(1),
        _shape{n0},
        _size(n0),
        _data((T*)malloc(sizeof(T)*_size)),
        _deleter([](T*p){free(p);})
        {};
    
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
    inline TArray<T>::TArray():_ndim(0),_deleter([](T*p){;}){};   
    

    template <typename T> 
    inline void TArray<T>::operator=(const T a) { for(index i=0;i<_size;i++) _data[i]=a;}

    template <typename T> 
    inline void TArray<T>::operator+=(const T a) { for(index i=0;i<_size;i++) _data[i]+=a;}

    template <typename T> 
    inline void TArray<T>::operator*=(const T a) { for(index i=0;i<_size;i++) _data[i]*=a;}

    template <typename T> 
    inline void TArray<T>::operator-=(const T a) { for(index i=0;i<_size;i++) _data[i]-=a;}

    template <typename T> 
    inline void TArray<T>::operator/=(const T a) { for(index i=0;i<_size;i++) _data[i]/=a;}

    template <typename T> 
    inline void TArray<T>::operator+=(const TArray<T> & a) { for(index i=0;i<_size;i++) _data[i]+=a._data[i];}

    template <typename T> 
    inline void TArray<T>::operator-=(const TArray<T> & a) { for(index i=0;i<_size;i++) _data[i]-=a._data[i];}

            
    template <typename T> 
    inline T TArray<T>::min() const { T min=_data[0]; for(index i=1;i<_size;i++) if (_data[i]<min) min=_data[i]; return min;}

    template <typename T> 
    inline T TArray<T>::max() const { T max=_data[0]; for(index i=1;i<_size;i++) if (_data[i]>max) max=_data[i]; return max;}

    template <typename T> 
    inline T TArray<T>::sum() const { T sum=_data[0]; for(index i=1;i<_size;i++) sum+=_data[i]; return sum;}

    template <typename T> 
    inline T TArray<T>::norm2() const  { return sqrt(dot(*this,*this));}

    template <typename T> 
    inline T TArray<T>::norm1() const  { T sum=std::abs(_data[0]); for(index i=1;i<_size;i++) sum+=std::abs(_data[i]); return sum; }

    template <typename T> 
    inline T TArray<T>::normi() const  { T x,max=std::abs(_data[0]); for(index i=1;i<_size;i++) if ((x=std::abs(_data[i]))>max) max=x; return max; }

    template <typename T> 
    inline void TArray<T>::fill(const TArray<T> &A) { for(index i=0;i<_size;i++) _data[i]=A._data[i];}

    template <typename T> 
    inline void TArray<T>::fill(std::function< T(const T)> f,const TArray<T> & A) { for(index i=0;i<A._size;i++) _data[i]=f(A._data[i]);}

    template <typename T> 
    inline void TArray<T>::operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B)  { for(index i=0;i<A._size;i++) f(A._data[i],B._data[i]);}

    template <typename T> 
    inline void TArray<T>::operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C)  { for(index i=0;i<A._size;i++) f(A._data[i],B._data[i],C._data[i]);}

    template <typename T> 
    inline T TArray<T>::dot(const TArray<T>& A,const TArray<T> &B) { T xdot=A._data[0]*B._data[0];for(index i=1;i<A._size;i++) xdot+=A._data[i]*B._data[i]; return xdot;}


}
