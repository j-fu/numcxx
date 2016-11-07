#ifndef NUMCXX_HXX
#define NUMCXX_HXX

namespace numcxx
{
    using index= unsigned int;
}

#include "expression.hxx"
#include "tarray.hxx"
#include "tarray1.hxx"
#include "tarray2.hxx"
#include "tmatrix.hxx"
#include "tsolverlapacklu.hxx"

/// Namespace for all parts of numcxx library
namespace numcxx
{
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using DMatrix=TMatrix<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
}
#endif

