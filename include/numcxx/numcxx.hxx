#ifndef NUMCXX_HXX
#define NUMCXX_HXX

/// Numcxx template library.
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
#include "tprecon-jacobi.hxx"
#include "util.hxx"

namespace numcxx
{
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using DMatrix=TMatrix<double>;
    using DSparseMatrix=TSparseMatrix<double>;
    using DSolverLapackLU=TSolverLapackLU<double>;
    using DSolverUMFPACK=TSolverUMFPACK<double>;
    using DPreconJacobi=TPreconJacobi<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;

  double norm1(const std::shared_ptr<DArray1> a) {return norm1(*a);}
    double norm2(const std::shared_ptr<DArray1> a) {return norm2(*a);}
    double normi(const std::shared_ptr<DArray1> a){return normi(*a);}

}


#endif

