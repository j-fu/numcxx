#ifndef NUMCXX_UTIL_HXX
#define NUMCXX_UTIL_HXX
/// untility functions

namespace  numcxx
{

    /// Maximum norm of array or expression
    ///
    /// \return Maximum ( \f$l^\infty\f$) norm.
    template <typename A> inline typename A::value_type normi(const A& a);

    /// Sum norm of array or expression.
    ///
    /// \return Sum ( \f$l^1\f$) norm.
    template <typename A> inline typename A::value_type norm1(const A& a);

    /// Euklidean norm of array or expression.
    ///
    /// \return Euklidean ( \f$l^2\f$) norm 
    template <typename A> inline typename A::value_type norm2(const A& a);

    /// Euklidean norm of array or expression.
    ///
    /// \return Euklidean ( \f$l^2\f$) norm 
    template <typename A> inline typename A::value_type norm(const A& a) { return norm2(a);}


    /// Dot product of array or expression.
    ///
    /// \param A First argument.
    /// \param B Second argument.
    /// \return  Dot product
    template <typename A, typename B> inline typename A::value_type dot(const A& a, const B&b);

    
    /// Minimum of of array or expression.
    ///
    /// \return Minimum of all elements in array.
    template <typename A> inline typename A::value_type min(const A&a);
    
    /// Maximum of  array or expression.
    ///
    /// \return Maximum of all elements in array.
    template <typename A> inline  typename A::value_type max(const A&a);
    
    /// Sum of array or expression.
    ///
    /// \return Sum of all elements in array.
    template <typename A> inline  typename A::value_type sum(const A&a);


    /// Create array of $n$ equally spaced  entries.
    template <typename T>
    inline std::shared_ptr<TArray1<T>>  linspace(T x0, T x1, int n);


    /// cpu time in seconds
    inline double cpu_clock();
    
    /// wall clock time in seconds
    inline double wall_clock();

      using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    inline double norm1(const std::shared_ptr<DArray1> a) {return norm1(*a);}
    inline double norm2(const std::shared_ptr<DArray1> a) {return norm2(*a);}
    inline double normi(const std::shared_ptr<DArray1> a){return normi(*a);}

}

#include "util.ixx"
#endif
