#ifndef NUMCXX_HXX
#define NUMCXX_HXX

#include "tarray.hxx"
#include "tarray1.hxx"
#include "tarray2.hxx"
#include "tmatrix.hxx"
#include "tsolverlapacklu.hxx"
#include "expression.hxx"

namespace numcxx
{
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using DMatrix=TMatrix<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
}
#endif

