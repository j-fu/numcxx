#include <ctime>

namespace  numcxx
{


  template <typename T>
  inline std::shared_ptr<TArray1<T>>  linspace(T x0, T x1, int n)
  { 
    auto pX=std::make_shared<TArray1<T>>(n);
    auto &X=*pX;
    T h=1.0/(T)(n-1);
    T x=0.0;
    for (int i=0;i<n;i++,x+=h) X[i]=x;
    return pX;
  }

  template <typename A> inline TArray1<typename A::value_type> arrayexpr(const A& a)
  {
    TArray1<typename A::value_type> v=a;
    return std::move(v);
  }


  template <typename A> inline typename A::value_type normi(const A& a)
  {
    typename A::value_type norm=std::abs(a[0]);
    for(index i=1; i<a.size(); i++ )
    {
      typename A::value_type x=std::abs(a[i]);
      if (x>norm) norm=x;
    }
    return norm;
  }

  template <typename A> inline typename A::value_type norm1(const A& a)
  {
    typename A::value_type norm=std::abs(a[0]);
    for(index i=1; i<a.size(); i++ )
    {
      norm+=std::abs(a[i]);
    }
    return norm;
  }

  template <typename A> inline typename A::value_type norm2(const A& a)
  {
    typename A::value_type norm=0.0;
    for(index i=0; i<a.size(); i++ )
    {
      typename A::value_type x=a[i];
      norm+=x*x;
    }
    return sqrt(norm);
  }


  template <typename A> inline typename A::value_type min(const A&a)
  {
    typename A::value_type m=a[0];
    for(index i=1; i<a.size(); i++ )
    {
      typename A::value_type x=a[i];
      if (x<m) m=x;
    }
    return m;
  }
    
  template <typename A> inline typename A::value_type max(const A&a)
  {
    typename A::value_type m=a[0];
    for(index i=1; i<a.size(); i++ )
    {
      typename A::value_type x=a[i];
      if (x>m) m=x;
    }
    return m;
  }

  template <typename A> inline  typename A::value_type sum(const A&a)
  {
    typename A::value_type s=a[0];
    for(index i=1; i<a.size(); i++ )
    {
      s+=a[i];
    }
    return s;
  }
    
  template <typename A, typename B> inline typename A::value_type dot(const A& a, const B&b)
  {
    typename A::value_type dot=0.0;
    for(index i=0; i<a.size(); i++ )
    {
      dot+=a[i]*b[i];
    }
    return dot;
  }    

  inline double cpu_clock() 
  {
    return (double)std::clock()/(double)CLOCKS_PER_SEC;
  }
    
  inline double wall_clock()
  { 
    time_t t;
    t=time(&t);
    return ((double)t);
  }

}
