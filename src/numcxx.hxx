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

extern "C"
{
  extern void dgetrf_(int *n, int *m, double *a, int *lda, int* ipiv, int *info);
  extern void dgetrs_(char *trans,int *n, const int *nrhs, double*a, int* lda, int *ipiv , double *b, int *ldb, int *info );
  extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int * N, int *K,  double *ALPHA, double *A, int *LDA,  double *B, int *LDB,  double *BETA,double * C, int *LDC);
}

namespace numcxx
{
    using index= unsigned int;


    template <class A> class array_ptr: public std::shared_ptr <A>
    {
        typedef typename A::base_type T;
        friend A;
    
    public:
        array_ptr(A *a):       std::shared_ptr<A>(a)    {};
        inline T & operator()(index i0) {return (*(std::shared_ptr<A>::get()))(i0);};
        inline T & operator()(index i0, index  i1) {return (*(std::shared_ptr<A>::get()))(i0,i1);};
        operator  A&()  const { return  *(std::shared_ptr<A>::get()); }

        array_ptr(){};
    };


    class TArrayBase 
    {
    public:
        TArrayBase() {}
        virtual ~TArrayBase() {}
    };
    
    template<typename T> class TArray1;  
    template<typename T> class TArray2;  

    template<typename T> class TArray: public TArrayBase
    {
        friend class TArray1<T>;
        friend class TArray2<T>;
    public:
        T*rawdata() { return data;}
        const index ndim;
        index size;
        std::vector<index> shape;

        T & operator()(index i0, index i1);
        T & operator()(index i0);

        void fill(T x);

        void operator=(const T a);
        void operator+=(const T a);
        void operator*=(const T a);
        void operator-=(const T a);
        void operator/=(const T a);
        T min() const { T min=data[0]; for(index i=1;i<size;i++) if (data[i]<min) min=data[i]; return min;}
        T max() const { T max=data[0]; for(index i=1;i<size;i++) if (data[i]>max) max=data[i]; return max;}
        T sum() const { T sum=data[0]; for(index i=1;i<size;i++) sum+=data[i]; return sum;}
        T norm2() const  { return dot(*this,*this);}
        T norm1() const  { T sum=std::abs(data[0]); for(index i=1;i<size;i++) sum+=abs(data[i]); return max; }
        T normi() const  { T x,max=std::abs(data[0]); for(index i=1;i<size;i++) if ((x=abs(data[i]))>max) max=x; return max; }
        void fill(const TArray<T> &A) { for(index i=0;i<size;i++) data[i]=A.data[i];}
        void fill(std::function< T(const T)> f,const TArray<T> & A) { for(index i=0;i<A.size;i++) data[i]=f(A.data[i]);}
        static void operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B)  { for(index i=0;i<A.size;i++) f(A.data[i],B.data[i]);}
        static void operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C)  { for(index i=0;i<A.size;i++) f(A.data[i],B.data[i],C.data[i]);}
        static T dot(const TArray<T>& A,const TArray<T> &B) { T xdot=A->data[0]*B->data[0];for(index i=1;i<A.size;i++) xdot+=A.data[i]*B.data[i]; return xdot;}
        typedef T base_type;

    private:
        const bool owndata=true;
        std::shared_ptr<TArrayBase>base =0;
        void check_bounds(index acc_dim, index acc_ndim, index acc_idx) const;

    protected:
        T* data=nullptr;
        index idx(index i0) const;
        index idx(index i0,index i1)  const;
        index idx(index i0,index i1,index i2)  const;
    
        TArray();
        ~TArray();

        TArray(index n0);
        TArray(index n0, T*aliendata);
        TArray(index n0, T*aliendata,std::shared_ptr<TArrayBase> xbase);

        TArray(index n0, index n1);
        TArray(index n0, index n1,T*aliendata);
        TArray(index n0, index n1, T*aliendata,std::shared_ptr<TArrayBase> xbase);
    };
  

    template<typename T> class TArray1: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;

        TArray1();
        TArray1(index n0);
        TArray1(index n0, T*aliendata);
        TArray1(index n0, T*aliendata, std::shared_ptr<TArrayBase> base);
        static array_ptr<TArray1 <T> > create(index n1);

        TArray1<T>(const std::initializer_list<T> &il ):TArray1(il.size())
        {
            index i=0;
            for (auto x = il.begin() ; x != il.end(); x++,i++) data[i]= *x;
        }
        static array_ptr<TArray1 <T> > create(const std::initializer_list<T> il)
        {
            auto A=new TArray1<T>(il.size());
            index i=0;
            for (auto x = il.begin() ; x != il.end(); x++,i++) A->data[i]= *x;
            return array_ptr<TArray1<T>> (A);
        }

        array_ptr<TArray1 <T> > copy() const;
        array_ptr<TArray1 <T> > clone() const;

        T & operator[](index i0);

        T item(index i0);
        void itemset(index i0, T x);

        T __getitem__(index i) const;
        void __setitem__(index i,T d);

        friend std::ostream & operator << (std::ostream & s, const array_ptr<TArray1<T>>&A)
        {
            return operator<<(s,*A);
        }


        friend std::ostream & operator << (std::ostream & s, TArray1<T> &A)
        {
            for (index i=0;i<A.size;i++) s <<"[" << i << "]: " <<A(i) << std::endl << std::flush;
            return s;
        }
        
    private:
        using TArray<T>::data;
        using TArray<T>::idx;

    };


    template<typename T> class TArray2: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;
        
        TArray2();
        TArray2(index n0, index n1);
        TArray2(index n0, index n1, T*aliendata);
        TArray2(index n0,index n1, T*aliendata, std::shared_ptr<TArrayBase> base);
        static array_ptr<TArray2 <T> > create(index n0,index n1);
        array_ptr<TArray2 <T> > copy() const;
        array_ptr<TArray2 <T> > clone() const;


        TArray2<T>(const  std::initializer_list<std::initializer_list<T>> &il ):TArray2(il.size(), il.begin()->size())
        {
            index i=0;
            
            for (auto jl = il.begin() ; jl != il.end(); jl++,i++)
            {
                index j=0;
                for (auto x = jl->begin() ; x != jl->end(); x++,j++) 
                    data[idx(i,j)]= *x;
            }
        }
        static array_ptr<TArray2 <T> > create(const  std::initializer_list<std::initializer_list<T>> &il)
        {
            auto A=new TArray2<T>(il.size(),il.begin()->size());
            index i=0;
            for (auto jl = il.begin() ; jl != il.end(); jl++,i++)
            {
                index j=0;
                for (auto x = jl->begin() ; x != jl->end(); x++,j++) 
                    (*A)(i,j)= *x;
            }
            return array_ptr<TArray2<T>> (A);
        }
        
        T * operator[](index i0);
        T item(index i0,index i1);
        void itemset(index i0, index i1, T x);
        array_ptr<TArray1 <T> > const __getitem__(index i0);

        friend std::ostream & operator << (std::ostream & s, const array_ptr<TArray2<T>>&A)
        {
            return operator<<(s,*A);
        }


        friend std::ostream & operator << (std::ostream & s, TArray2<T> &A)
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
        
        static void lu_decomp(TArray2<T> &lu, TArray1<int> &ipiv) 
        {
            if (typeid(T)==typeid(double))
            {
                int n=lu.shape[0];
                int info;
                dgetrf_(&n,&n,lu.rawdata(),&n,ipiv.rawdata(),&info);
            }
        }

        std::tuple<array_ptr<TArray2<T>>, array_ptr<TArray1<int>>> lu_decomp() const
        {
            assert(shape[0]==shape[1]);
            auto lu=this->copy();
            auto ipiv=TArray1<int>::create(shape[0]);
            lu_decomp(*lu,*ipiv);
            return std::make_tuple(lu,ipiv);
        }
        static array_ptr<TArray1<T>> lu_solve(TArray2<T> &lu, TArray1<int> &ipiv,const TArray1<T> &rhs)
        {
            auto sol=rhs.clone();
            lu_solve(lu,ipiv,sol,rhs);
            return sol;
        }
        static void lu_solve(TArray2<T> &lu, TArray1<int> &ipiv,TArray1<T> &sol, const TArray1<T> &rhs)
        {
            sol.fill(rhs);
            if (typeid(T)==typeid(double))
            {
              char trans[2]={'T','\0'};
              int n=lu.shape[0];
              int one=1;
              int info;
              dgetrs_(trans,&n,&one,lu.rawdata(),&n,ipiv.rawdata(),sol.data,&n,&info);
            }            
        }
        void solve(TArray1<T> &sol, const TArray1<T> &rhs) const
        {
            array_ptr< TArray2<T> > lu;
            array_ptr< TArray1<int> > ipiv;
            std::tie(lu,ipiv)=lu_decomp();
            lu_solve(*lu,*ipiv,sol,rhs);
        }

        array_ptr<TArray1<T>> solve(const TArray1<T> &rhs) const
        {
            auto sol=rhs.clone();
            solve(sol,rhs);
            return sol;
        }


        void apply(const TArray1<T> &u, TArray1<T> &v)
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
                       data,&n,
                       u.data,&n,
                       &dzero,
                       v.data,&n);
            }
        }
        array_ptr<TArray1<T>> apply(const TArray1<T> &u)
        {
            auto v=u.clone();
            apply(u,*v);
            return v;
        }

        
    private:
        using TArray<T>::data;
        using TArray<T>::idx;
    };
    
    
    
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
    
/////////////////////////////////////////////////
// Inline methods of TArray<T>

    template <typename T> 
    inline T & TArray<T>::operator()(index i0) { return data[idx(i0)];};


    template <typename T> 
    inline T & TArray<T>::operator()(index i0, index i1) { return data[idx(i0,i1)];};

    template <typename T> 
    inline void TArray<T>::fill(T x)
    {
        for (index i=0;i<size;i++) data[i]=x;
    }
    
    
    template <typename T> 
    inline void TArray<T>::check_bounds(index acc_dim, index acc_ndim, index acc_idx) const
    {
        if (acc_ndim!=ndim) 
        {
            char errormsg[80];
            snprintf(errormsg,80,"numcxx::TArray::check_bounds: attempt of %uD access of %uD array",acc_ndim,ndim);
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
        ndim(1),
        shape{n0},
        size(n0)
        {
            data=new T[size];
        };
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, index n1):
        ndim(2),
        shape{n0,n1},
        size(n0*n1)
        {
            data=new T[size];
        };
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, T*aliendata):
        data(aliendata),
        ndim(1),
        shape{n0},
        size(n0),
        owndata(false){};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, T*aliendata,std::shared_ptr<TArrayBase> xbase):
        TArray(n0,aliendata){base=xbase;};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, index n1,T*aliendata):
        data(aliendata),
        ndim(2),
        shape{n0,n1},
        size(n0*n1),
        owndata(false){};
    
    template <typename T> 
    inline TArray<T>::TArray(index n0, index n1, T*aliendata,std::shared_ptr<TArrayBase> xbase):
        TArray(n0,n1,aliendata){base=xbase;};
    
    template <typename T> 
    inline TArray<T>::~TArray()
    {
        if (owndata) delete[] data;
    };
    
    template <typename T> 
    inline TArray<T>::TArray():ndim(0),owndata(false){};   
    

    template <typename T> 
    inline void TArray<T>::operator=(const T a) { for(index i=0;i<size;i++) data[i]=a;}
    template <typename T> 
    inline void TArray<T>::operator+=(const T a) { for(index i=0;i<size;i++) data[i]+=a;}
    template <typename T> 
    inline void TArray<T>::operator*=(const T a) { for(index i=0;i<size;i++) data[i]*=a;}
    template <typename T> 
    inline void TArray<T>::operator-=(const T a) { for(index i=0;i<size;i++) data[i]-=a;}
    template <typename T> 
    inline void TArray<T>::operator/=(const T a) { for(index i=0;i<size;i++) data[i]/=a;}


/////////////////////////////////////////////////
// Inline methods of TArray1<T>

    template <typename T> 
    inline TArray1<T>::TArray1():TArray<T>(){};   

    template <typename T> 
    inline TArray1<T>::TArray1(index n0):TArray<T>(n0){}

    template <typename T> 
    inline TArray1<T>::TArray1(index n0, T*aliendata):
        TArray<T>(n0,aliendata){}
    
    template <typename T> 
    inline TArray1<T>::TArray1(index n0, T*aliendata, std::shared_ptr<TArrayBase> base):
        TArray<T>(n0,aliendata,base){}; 
    
    template <typename T> 
    inline array_ptr<TArray1 <T> > TArray1<T>::create(index n1)
    {
        return array_ptr<TArray1 <T> >(new TArray1<T>(n1));
    }
    
    template <typename T> 
    inline T & TArray1<T>::operator[](index i0) { return data[idx(i0)];};
    
    template <typename T> 
    inline T TArray1<T>::item(index i0) { return data[idx(i0)];};

    template <typename T> 
    inline void TArray1<T>::itemset(index i0, T x) { data[idx(i0)]=x;};

    template <typename T> 
    inline array_ptr<TArray1 <T> > TArray1<T>::copy() const
    {
        auto x=create(shape[0]);
        for (index i=0;i<size;i++) x->data[i]=data[i];
        return x;
    }

    template <typename T> 
    inline array_ptr<TArray1 <T> > TArray1<T>::clone() const
    {
        return create(shape[0]);
    }

    template <typename T> 
    inline T TArray1<T>::__getitem__(index i) const { return data[idx(i)]; } 

    template <typename T> 
    inline void TArray1<T>::__setitem__(index i,T d) { data[idx(i)]=d; } 


/////////////////////////////////////////////////
// inline methods of TArray2<T>

    
    template <typename T> 
    inline TArray2<T>::TArray2():TArray<T>(){};   
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0, index n1):TArray<T>(n0,n1) {}
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0, index n1, T*aliendata):
        TArray<T>(n0,n1,aliendata){}
    
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0,index n1, T*aliendata, std::shared_ptr<TArrayBase> base):
        TArray<T>(n0,n1,aliendata,base) {}
    
    template <typename T> 
    inline array_ptr<TArray2 <T> > TArray2<T>::create(index n0,index n1)
    {
        return array_ptr<TArray2 <T> >(new TArray2(n0,n1));
    }
    
    
    template <typename T> 
    inline T * TArray2<T>::operator[](index i0) { return &data[idx(i0,0)];};

    template <typename T> 
    inline T TArray2<T>::item(index i0,index i1) { return data[idx(i0,i1)];};

    template <typename T> 
    inline void TArray2<T>::itemset(index i0, index i1, T x) { data[idx(i0,i1)]=x;};

    template <typename T> 
    inline array_ptr<TArray2 <T> > TArray2<T>::copy() const
    {
        auto x=create(shape[0],shape[1]);
        for (index i=0;i<size;i++) x->data[i]=data[i];
        return x;
    }

    template <typename T> 
    inline array_ptr<TArray2 <T> > TArray2<T>::clone() const
    {
        return create(shape[0],shape[1]);
    }

    template <typename T> 
    inline array_ptr<TArray1 <T> > const TArray2<T>::__getitem__(index i0) 
    { 
        return array_ptr<TArray1<T>>(new TArray1<T>(shape[1], &data[idx(i0,0)]));
    } 
    

}

#endif

