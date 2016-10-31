/* A header-only, numpy compatible, lightweight multidimensional array class for C++11.
   
   In  difference to  ndarray (https://github.com/ndarray/ndarray)  it
   does not  need boost, but  relies on swig  for the creation  of the
   python binding.   It provides reference-counted  conversion from/to
   numpy  arrays  without  copying  and  array_ptr,  a  smart  pointer
   derieved from std::shared_ptr.
   
   Method names are inspired by numpy, but by far not comprehensive.
   
   row major layout!
   
   http://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays/  
   
   pybind11? Young project.
   https://community.lsst.org/t/using-pybind11-instead-of-swig-to-wrap-c-code/1096

*/
#ifndef NUMCXX_HXX
#define NUMCXX_HXX

#include <vector>
#include <typeinfo>
#include <memory>
#include <stdexcept> 
#include <ostream>
#include <tuple>

#include <cmath>
#include <cassert>
#include <cstdlib>


namespace numcxx
{
    using index= unsigned int;
    
    template<typename T> class TArray
    {
    public:
        T*data() const { return _data;}
        std::vector<index> shape;
        index ndim() {return _ndim;}
        index size() {return _size;}

        T & operator()(index i0, index i1);
        T & operator()(index i0);

        void fill(T x);

        void operator=(const T a);
        void operator+=(const T a);
        void operator*=(const T a);
        void operator-=(const T a);
        void operator/=(const T a);
        T min() const;
        T max() const;
        T sum() const;
        T norm2() const;
        T norm1() const;
        T normi() const;
        void fill(const TArray<T> &A);
        void fill(std::function< T(const T)> f,const TArray<T> & A);
        static void operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B);
        static void operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C);
        static T dot(const TArray<T>& A,const TArray<T> &B);

        typedef T value_type;

    private:
        std::shared_ptr<void>base =0;
        void check_bounds(index acc_dim, index acc_ndim, index acc_idx) const;
        const std::function<void(T*p)> deleter;

    protected:
        const index _ndim;
        index _size;
        T* _data=nullptr;

        index idx(index i0) const;
        index idx(index i0,index i1)  const;
        index idx(index i0,index i1,index i2)  const;
    

        TArray(index n0);
        TArray(index n0, T*alien_data, std::function<void(T*p)> deleter);
        TArray(index n0, T*alien_data,std::shared_ptr<void> xbase);

        TArray(index n0, index n1);
        TArray(index n0, index n1,T*alien_data,std::function<void(T*p)> deleter);
        TArray(index n0, index n1, T*alien_data,std::shared_ptr<void> xbase);

        TArray();
        ~TArray();

    };


    template<typename T>
    class TArray1;


    template<typename T>
    inline std::ostream & operator << (std::ostream & s, TArray1<T> &A);

    template<typename T> class TArray1: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;

        TArray1();
        TArray1(index n0);
        TArray1(index n0, T*alien_data,std::function<void(T*p)> deleter);
        TArray1(index n0, T*alien_data, std::shared_ptr<void> base);
        TArray1(std::shared_ptr<std::vector<T>> v);
        TArray1(const std::initializer_list<T> &il );
        static std::shared_ptr<TArray1 <T> > create(index n1);
        static std::shared_ptr<TArray1 <T> > create(const std::initializer_list<T> il);

        std::shared_ptr<TArray1 <T> > copy() const;
        std::shared_ptr<TArray1 <T> > clone() const;

        T & operator[](index i0);

        T item(index i0);
        void itemset(index i0, T x);

        T __getitem__(index i) const;
        void __setitem__(index i,T d);

        friend std::ostream & operator<< <T>(std::ostream & s, TArray1<T> &A);

        
    private:
        using TArray<T>::_data;
        using TArray<T>::idx;
        using TArray<T>::_size;

    };



    template<typename T>
    class TArray2;

    template<typename T>
    inline std::ostream & operator << (std::ostream & s, TArray2<T> &A);

    template<typename T> 
    class TArray2: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;
        
        TArray2();
        TArray2(index n0, index n1);
        TArray2(index n0, index n1, T*alien_data, std::function<void(T*p)> deleter);
        TArray2(index n0,index n1, T*alien_data, std::shared_ptr<void> base);
        TArray2(const  std::initializer_list<std::initializer_list<T>> &il );
        static std::shared_ptr<TArray2 <T> > create(index n0,index n1);
        static std::shared_ptr<TArray2 <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);
        std::shared_ptr<TArray2 <T> > copy() const;
        std::shared_ptr<TArray2 <T> > clone() const;

        T * operator[](index i0);
        T item(index i0,index i1);
        void itemset(index i0, index i1, T x);
        std::shared_ptr<TArray1 <T> > const __getitem__(index i0);


        friend std::ostream & operator<< <T>(std::ostream & s, TArray2<T> &A);

    protected:
        using TArray<T>::_data;
        using TArray<T>::idx;
        using TArray<T>::_size;
    };
    

    template<typename T> 
    class TMatrix: public TArray2<T>
    {
    public:
        using TArray2<T>::size;
        using TArray2<T>::shape;

        TMatrix();
        TMatrix(index n0, index n1);
        TMatrix(index n0, index n1, T*alien_data, std::function<void(T*p)> deleter);
        TMatrix(index n0,index n1, T*alien_data, std::shared_ptr<void> base);
        TMatrix(const  std::initializer_list<std::initializer_list<T>> &il );
        static std::shared_ptr<TMatrix <T> > create(index n0,index n1);
        static std::shared_ptr<TMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);
        std::shared_ptr<TMatrix <T> > copy() const;
        std::shared_ptr<TMatrix <T> > clone() const;

        static void lu_decomp(TMatrix<T> &lu, TArray1<int> &ipiv); 
        static void lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,TArray1<T> &sol, const TArray1<T> &rhs);
        void solve(TArray1<T> &sol, const TArray1<T> &rhs) const;
        void apply(const TArray1<T> &u, TArray1<T> &v);

        static void inverse(TMatrix<T> &lu, TArray1<int> &ipiv, TMatrix<T>& inverse);
        void inverse(TMatrix<T>& inverse) const;

        std::tuple<std::shared_ptr<TMatrix<T>>, std::shared_ptr<TArray1<int>>> lu_decomp() const;
        static std::shared_ptr<TArray1<T>> lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,const TArray1<T> &rhs);
        std::shared_ptr<TArray1<T>> solve(const TArray1<T> &rhs) const;
        std::shared_ptr<TArray1<T>> apply(const TArray1<T> &u);

        static std::shared_ptr<TMatrix<T>> inverse(const TMatrix<T> &lu, const TArray1<int> &ipiv);
        std::shared_ptr<TMatrix<T>> inverse() const;
        
        


    private:
        using TArray2<T>::_data;
        using TArray2<T>::_size;
        using TArray2<T>::idx;

    };
    
    


    using DMatrix=TMatrix<double>;
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
    
/////////////////////////////////////////////////
// Inline methods of TArray<T>

    template <typename T> 
    inline T & TArray<T>::operator()(index i0) { return _data[idx(i0)];};


    template <typename T> 
    inline T & TArray<T>::operator()(index i0, index i1) { return _data[idx(i0,i1)];};

    template <typename T> 
    inline void TArray<T>::fill(T x)
    {
        for (index i=0;i<_size;i++) _data[i]=x;
    }
    
    
    template <typename T> 
    inline void TArray<T>::check_bounds(index acc_dim, index acc_ndim, index acc_idx) const
    {
        if (acc_ndim!=_ndim) 
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TArray::check_bounds: attempt of %uD access of %uD array",acc_ndim,_ndim);
            throw std::out_of_range(errormsg);
        }
        if ((acc_idx<0) || (acc_idx>=shape[acc_dim]))
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TArray::check_bounds: shape[%u]=%u but i%u=%u",acc_dim,shape[acc_dim],acc_dim,acc_idx);
            throw std::out_of_range(errormsg);
        }
    }  
    
    template <typename T> 
    inline index TArray<T>::idx(index i0) const
    {
#ifdef NUMCXX_CHECK_BOUNDS
        check_bounds(0,1,i0);
#endif
        return i0;
    }

    template <typename T> 
    inline index TArray<T>::idx(index i0,index i1)  const
    { 
#ifdef NUMCXX_CHECK_BOUNDS
        check_bounds(0,2,i0);
        check_bounds(1,2,i1);
#endif
        return i0*shape[1]+i1;
    }
    
    template <typename T> 
    inline index TArray<T>::idx(index i0,index i1, index i2)  const
    { 
#ifdef NUMCXX_CHECK_BOUNDS
        check_bounds(0,3,i0);
        check_bounds(1,3,i1);
        check_bounds(2,3,i2);
#endif
         return (i0*shape[0]+i1)*shape[1]+i2;
    }
    
    template <typename T> 
        inline TArray<T>::TArray(index n0):
        _ndim(1),
        shape{n0},
        _size(n0),
        _data((T*)malloc(sizeof(T)*_size)),
        deleter([](T*p){free(p);})
        {};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, index n1):
        _ndim(2),
        shape{n0,n1},
        _size(n0*n1),
        _data((T*)malloc(sizeof(T)*_size)),
        deleter([](T*p){free(p);})
        {};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, T*alien_data,std::function<void(T*p)> deleter):
        _data(alien_data),
        _ndim(1),
        shape{n0},
        _size(n0),
        deleter(deleter)
        {};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, T*alien_data,std::shared_ptr<void> xbase):
        TArray(n0,alien_data,[](T*p){;}){base=xbase;};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, index n1,T*alien_data, std::function<void(T*p)> deleter):
        _data(alien_data),
        _ndim(2),
        shape{n0,n1},
        _size(n0*n1),
        deleter(deleter)
        {};



    
    template <typename T> 
    inline TArray<T>::TArray(index n0, index n1, T*alien_data,std::shared_ptr<void> xbase):
        TArray(n0,n1,alien_data, [](T*p){;})
    {base=xbase;};
    
    template <typename T> 
    inline TArray<T>::~TArray()
    {
       deleter(_data);
    };
    
    template <typename T> 
    inline TArray<T>::TArray():_ndim(0),deleter([](T*p){;}){};   
    

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
    inline T TArray<T>::min() const { T min=_data[0]; for(index i=0;i<_size;i++) if (_data[i]<min) min=_data[i]; return min;}

    template <typename T> 
    inline T TArray<T>::max() const { T max=_data[0]; for(index i=0;i<_size;i++) if (_data[i]>max) max=_data[i]; return max;}

    template <typename T> 
    inline T TArray<T>::sum() const { T sum=_data[0]; for(index i=0;i<_size;i++) sum+=_data[i]; return sum;}

    template <typename T> 
    inline T TArray<T>::norm2() const  { return dot(*this,*this);}

    template <typename T> 
    inline T TArray<T>::norm1() const  { T sum=std::abs(_data[0]); for(index i=0;i<_size;i++) sum+=abs(_data[i]); return max; }

    template <typename T> 
    inline T TArray<T>::normi() const  { T x,max=std::abs(_data[0]); for(index i=0;i<_size;i++) if ((x=abs(_data[i]))>max) max=x; return max; }

    template <typename T> 
    inline void TArray<T>::fill(const TArray<T> &A) { for(index i=0;i<_size;i++) _data[i]=A._data[i];}

    template <typename T> 
    inline void TArray<T>::fill(std::function< T(const T)> f,const TArray<T> & A) { for(index i=0;i<A._size;i++) _data[i]=f(A._data[i]);}

    template <typename T> 
    inline void TArray<T>::operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B)  { for(index i=0;i<A._size;i++) f(A._data[i],B._data[i]);}

    template <typename T> 
    inline void TArray<T>::operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C)  { for(index i=0;i<A._size;i++) f(A._data[i],B._data[i],C._data[i]);}

    template <typename T> 
    inline T TArray<T>::dot(const TArray<T>& A,const TArray<T> &B) { T xdot=A->_data[0]*B->_data[0];for(index i=0;i<A._size;i++) xdot+=A._data[i]*B._data[i]; return xdot;}

/////////////////////////////////////////////////
// Inline methods of TArray1<T>

    template <typename T> 
    inline TArray1<T>::TArray1():TArray<T>(){};   

    template <typename T> 
    inline TArray1<T>::TArray1(index n0):TArray<T>(n0){}

    template <typename T> 
    inline TArray1<T>::TArray1(index n0, T*alien_data,std::function<void(T*p)> deleter):
        TArray<T>(n0,alien_data,deleter){}
    
    template <typename T> 
    inline TArray1<T>::TArray1(index n0, T*alien_data, std::shared_ptr<void> base):
        TArray<T>(n0,alien_data,base){}; 

    template <typename T> 
    inline TArray1<T>::TArray1(std::shared_ptr<std::vector<T>> v):
        TArray<T>(v->size(),v->data(),v){};


    

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::create(index n1)
    {
        return std::make_shared<TArray1 <T> >(n1);
    }
    
    template<typename T> 
    inline TArray1<T>::TArray1(const std::initializer_list<T> &il ):TArray1(il.size())
    {
        index i=0;
        for (auto x = il.begin() ; x != il.end(); x++,i++) _data[i]= *x;
    }
    
    template<typename T>
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::create(const std::initializer_list<T> il)
    {
        return std::make_shared<TArray1 <T> >(il);
    }
    
    template <typename T> 
    inline T & TArray1<T>::operator[](index i0) { return _data[idx(i0)];};
    
    template <typename T> 
    inline T TArray1<T>::item(index i0) { return _data[idx(i0)];};

    template <typename T> 
    inline void TArray1<T>::itemset(index i0, T x) { _data[idx(i0)]=x;};

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::copy() const
    {
        auto x=create(shape[0]);
        for (index i=0;i<_size;i++) x->_data[i]=_data[i];
        return x;
    }

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::clone() const
    {
        return create(shape[0]);
    }

    template <typename T> 
    inline T TArray1<T>::__getitem__(index i) const { return _data[idx(i)]; } 

    template <typename T> 
    inline void TArray1<T>::__setitem__(index i,T d) { _data[idx(i)]=d; } 

    template<typename T> 
    inline std::ostream & operator << (std::ostream & s, const std::shared_ptr<TArray1<T>>&A)
    {
        return operator<<(s,*A);
    }
    
    
    template<typename T> 
    inline std::ostream & operator << (std::ostream & s, TArray1<T> &A)
    {
        for (index i=0;i<A._size;i++) s <<"[" << i << "]: " <<A(i) << std::endl << std::flush;
        return s;
    }


/////////////////////////////////////////////////
// inline methods of TArray2<T>

    
    template <typename T> 
    inline TArray2<T>::TArray2():TArray<T>(){};   
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0, index n1):TArray<T>(n0,n1) {}
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0, index n1, T*alien_data,std::function<void(T*p)> deleter):
        TArray<T>(n0,n1,alien_data,deleter(deleter)){}
    
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0,index n1, T*alien_data, std::shared_ptr<void> base):
        TArray<T>(n0,n1,alien_data,base) {}
    
    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::create(index n0,index n1)
    {
        return std::make_shared<TArray2 <T> >(n0,n1);
    }
    
    
    template <typename T> 
    inline T * TArray2<T>::operator[](index i0) { return &_data[idx(i0,0)];};

    template <typename T> 
    inline T TArray2<T>::item(index i0,index i1) { return _data[idx(i0,i1)];};

    template <typename T> 
    inline void TArray2<T>::itemset(index i0, index i1, T x) { _data[idx(i0,i1)]=x;};

    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::copy() const
    {
        auto x=create(shape[0],shape[1]);
        for (index i=0;i<_size;i++) x->_data[i]=_data[i];
        return x;
    }

    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::clone() const
    {
        return create(shape[0],shape[1]);
    }

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > const TArray2<T>::__getitem__(index i0) 
    { 
        return std::shared_ptr<TArray1<T>>(new TArray1<T>(shape[1], &_data[idx(i0,0)], [](T*p){;}));
    } 



    template <typename T> 
    inline TArray2<T>::TArray2(const  std::initializer_list<std::initializer_list<T>> &il ):TArray2(il.size(), il.begin()->size())
    {
        index i=0;
        
        for (auto jl = il.begin() ; jl != il.end(); jl++,i++)
        {
            index j=0;
            for (auto x = jl->begin() ; x != jl->end(); x++,j++) 
                _data[idx(i,j)]= *x;
        }
    }
    

    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::create(const  std::initializer_list<std::initializer_list<T>> &il)
    {
        return std::make_shared<TArray2<T>> (il);
    }
        

    
    
    template <typename T> 
    inline std::ostream & operator << (std::ostream & s, TArray2<T> &A)
    {
        s << "    ";
        for (index j=0;j<A.shape[1];j++) 
            s << "[" << j << "]     ";
        s<< std::endl;
        for (index i=0;i<A.shape[0];i++) 
        {
            s << "[" << i << "]: ";
            for (index j=0;j<A.shape[1];j++) 
                s << A(i,j) << "   ";
            s<< std::endl;
            }
        s << std::flush;
        return s;
    }
    
    template <typename T> 
    inline std::ostream & operator << (std::ostream & s, const std::shared_ptr<TArray2<T>>&A)
    {
        return operator<<(s,*A);
    }

///////////////////////////////////////////////
/// inline methods of TMatrix
    extern "C"
    {
        extern void dgetrf_(int *n, int *m, double *a, int *lda, int* ipiv, int *info);
        extern void dgetrs_(char *trans,int *n, const int *nrhs, double*a, int* lda, int *ipiv , double *b, int *ldb, int *info );
        extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  double *ALPHA, const double *A, int *LDA,  double *B, int *LDB,  double *BETA,double * C, int *LDC);
    }


    template <typename T> 
    inline TMatrix<T>::TMatrix(): TArray2<T>(){};

    template <typename T> 
    inline TMatrix<T>::TMatrix(index n0, index n1): TArray2<T>(n0,n1){};

    template <typename T> 
    inline TMatrix<T>::TMatrix(index n0, index n1, T*alien_data,std::function<void(T*p)> deleter): 
        TArray2<T>(n0,n1,alien_data,deleter){};

    template <typename T> 
    inline TMatrix<T>::TMatrix(index n0,index n1, T*alien_data, std::shared_ptr<void> base): 
        TArray2<T>(n0,n1,alien_data,base){};

    template <typename T> 
    inline TMatrix<T>::TMatrix(const  std::initializer_list<std::initializer_list<T>> &il ): TArray2<T>(il){};

    template <typename T> 
    inline std::shared_ptr<TMatrix <T> > TMatrix<T>::create(index n0,index n1)
    {
        return std::make_shared<TMatrix<T>> (n0,n1);
    }

    template <typename T> 
    inline  std::shared_ptr<TMatrix <T> > TMatrix<T>::create(const  std::initializer_list<std::initializer_list<T>> &il)
    {
        return std::make_shared<TMatrix<T>> (il);
    }

    template <typename T> 
    inline  std::shared_ptr<TMatrix <T> > TMatrix<T>::copy() const
    {
        auto x=create(shape[0],shape[1]);
        for (index i=0;i<_size;i++) x->_data[i]=_data[i];
        return x;
    }

    template <typename T> 
    inline   std::shared_ptr<TMatrix <T> > TMatrix<T>::clone() const
    { 
        return create(shape[0],shape[1]);
    } 

    
    
    template <typename T> 
    inline   void TMatrix<T>::lu_decomp(TMatrix<T> &lu, TArray1<int> &ipiv) 
    {
        if (typeid(T)==typeid(double))
        {
            int n=lu.shape[0];
            int info;
            dgetrf_(&n,&n,lu.data(),&n,ipiv.data(),&info);
        }
    }
    
    template <typename T> 
    inline   std::tuple<std::shared_ptr<TMatrix<T>>, std::shared_ptr<TArray1<int>>> TMatrix<T>::lu_decomp() const
    {
        assert(shape[0]==shape[1]);
        auto lu=this->copy();
        auto ipiv=TArray1<int>::create(shape[0]);
        lu_decomp(*lu,*ipiv);
        return std::make_tuple(lu,ipiv);
    }

    template <typename T> 
    inline   std::shared_ptr<TArray1<T>> TMatrix<T>::lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,const TArray1<T> &rhs)
    {
        auto sol=rhs.clone();
        lu_solve(lu,ipiv,*sol,rhs);
        return sol;
    }

    template <typename T> 
    inline   void TMatrix<T>::lu_solve(TMatrix<T> &lu, TArray1<int> &ipiv,TArray1<T> &sol, const TArray1<T> &rhs)
    {
        sol.fill(rhs);
        if (typeid(T)==typeid(double))
        {
            char trans[2]={'T','\0'};
            int n=lu.shape[0];
            int one=1;
            int info;
            dgetrs_(trans,&n,&one,lu.data(),&n,ipiv.data(),sol.data(),&n,&info);
            }            
        }

    template <typename T> 
    inline   void TMatrix<T>::solve(TArray1<T> &sol, const TArray1<T> &rhs) const
        {
            std::shared_ptr< TMatrix<T> > lu;
            std::shared_ptr< TArray1<int> > ipiv;
            std::tie(lu,ipiv)=lu_decomp();
            lu_solve(*lu,*ipiv,sol,rhs);
        }
    
    template <typename T> 
    inline   std::shared_ptr<TArray1<T>> TMatrix<T>::solve(const TArray1<T> &rhs) const
        {
            auto sol=rhs.clone();
            solve(*sol,rhs);
            return sol;
        }
    

    template <typename T> 
    inline   void TMatrix<T>::apply(const TArray1<T> &u, TArray1<T> &v)
        {
            if (typeid(T)==typeid(double))
            {
                char transmat[2]={'T','\0'};
                char transvec[2]={'N','\0'};
                int n=shape[0];
                int ione=1;
                double done=1.0;
                double dzero=0.0;
                
                dgemm_(transmat,
                       transvec,
                       &n,&ione,&n,
                       &done,
                       _data,&n,
                       u.data(),&n,
                       &dzero,
                       v.data(),&n);
            }
        }

    template <typename T> 
    inline   std::shared_ptr<TArray1<T>> TMatrix<T>::apply(const TArray1<T> &u)
    {
        auto v=u.clone();
        apply(u,*v);
        return v;
    }
    

}

#endif

