/* 
   A header-only, lightweight multidimensional array class for C++11.
   
   In  difference to  ndarray (https://github.com/ndarray/ndarray)  it
   does not  need boost, but  relies on swig  for the creation  of the
   python binding.   It provides reference-counted  conversion from/to
   numpy  arrays  without  copying.
   
   Method names are inspired by numpy, but by far not comprehensive.
*/

#ifndef NUMCXX_HXX
#define NUMCXX_HXX

#include "tarray.hxx"
#include "tarray1.hxx"
#include "tarray2.hxx"
#include "tmatrix.hxx"

namespace numcxx
{
    using DMatrix=TMatrix<double>;
    using DArray1=TArray1<double>;
    using DArray2=TArray2<double>;
    using IArray1=TArray1<int>;
    using IArray2=TArray2<int>;
}
#endif

