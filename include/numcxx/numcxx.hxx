#ifndef NUMCXX_HXX
#define NUMCXX_HXX

#include "tarray.hxx"
#include "tarray1.hxx"
#include "tarray2.hxx"
#include "tmatrix.hxx"
#include "expression.hxx"

namespace numcxx
{
    using DMatrix=TMatrix<double>;
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
}
#endif

