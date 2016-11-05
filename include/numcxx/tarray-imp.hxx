/////////////////////////////////////////////////
// Inline methods of TArray<T>
#include <cstdlib>
#include <cmath>

namespace numcxx
{


    template <typename T>
    inline const T& TArray<T>::operator[](const index i) const
    {
        return _data[i];
    }

    template <typename T, typename EXPR,
              typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
    inline void assign(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
    {
//            resize( expr.size() );
        T *data=A.data();
        for(index i=0; i<expr.size(); i++ ) data[i] = expr[i];
    }

    template <typename T, typename VAL,
              typename= typename std::enable_if<!std::is_class<VAL>::value, VAL>::type>
    inline void assign(TArray<T>& A, const  VAL& a) 
    {
        T *data=A.data();
        for(index i=0; i<A.size(); i++ ) data[i] = a;
    }

    template <typename T, typename EXPR,
              typename= typename std::enable_if<std::is_class<EXPR>::value, EXPR>::type>
    inline void xadd(TArray<T>& A, const  EXPR& expr , const EXPR *x=0) 
    {
//            resize( expr.size() );
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
//            resize( expr.size() );
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
//            resize( expr.size() );
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
//            resize( expr.size() );
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


    template <typename A> double normi(const A& a)
    {
        double norm=std::abs(a[0]);
        for(index i=1; i<a.size(); i++ )
        {
            double x=std::abs(a[i]);
            if (x>norm) norm=x;
        }
        return norm;
    }

    template <typename A> double norm1(const A& a)
    {
        double norm=std::abs(a[0]);
        for(index i=1; i<a.size(); i++ )
        {
           norm+=std::abs(a[i]);
        }
        return norm;
    }

    template <typename A> double norm2(const A& a)
    {
        double norm=0.0;
        for(index i=0; i<a.size(); i++ )
        {
            double x=a[i];
            norm+=x*x;
        }
        return sqrt(norm);
    }


    template <typename A> double min(const A&a)
    {
        double m=a[0];
        for(index i=1; i<a.size(); i++ )
        {
            double x=a[i];
            if (x<m) m=x;
        }
        return m;
    }

    template <typename A> double max(const A&a)
    {
        double m=a[0];
        for(index i=1; i<a.size(); i++ )
        {
            double x=a[i];
            if (x>m) m=x;
        }
        return m;
    }

    template <typename A>  double sum(const A&a)
    {
        double s=a[0];
        for(index i=1; i<a.size(); i++ )
        {
           s+=a[i];
        }
        return s;
    }

    template <typename A, typename B> double dot(const A& a, const B&b)
    {
        double dot=0.0;
        for(index i=0; i<a.size(); i++ )
        {
            dot+=a[i]*b[i];
        }
        return dot;
    }

    template <typename T> 
    inline        T*TArray<T>::data() const { return _data;}

    template <typename T> 
    inline    index TArray<T>::ndim() const {return _ndim;}

    template <typename T> 
    inline    index TArray<T>::size() const {return _size;}

    template <typename T> 
    inline    index TArray<T>::shape(const index dim)  const {return _shape[dim];}
    
    
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
    inline void TArray<T>::operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B)  { for(index i=0;i<A._size;i++) f(A._data[i],B._data[i]);}

    template <typename T> 
    inline void TArray<T>::operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C)  { for(index i=0;i<A._size;i++) f(A._data[i],B._data[i],C._data[i]);}


}
