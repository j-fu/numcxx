#ifndef NUMCXX_HXX
#define NUMCXX_HXX

namespace numcxx
{
    using index= unsigned int;
}

#include "expression.ixx"
#include "tarray.hxx"
#include "tarray1.hxx"
#include "tarray2.hxx"
#include "tmatrix.hxx"
#include "tsparsematrix.hxx"
#include "tsolver-lapacklu.hxx"
#include "tsolver-umfpack.hxx"
#include "util.hxx"

/// Namespace for all parts of numcxx library
namespace numcxx
{
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using DMatrix=TMatrix<double>;
    using DSparseMatrix=TSparseMatrix<double>;
    using DSolverLapackLU=TSolverLapackLU<double>;
    using DSolverUMFPACK=TSolverUMFPACK<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
}
#endif

