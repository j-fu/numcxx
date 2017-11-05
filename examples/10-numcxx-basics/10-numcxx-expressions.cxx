///
/// \example 10-numcxx-expressions.cxx
///
/// Demonstrate the basic usage of vectors and expression templates in numcxx.
/// 

#include <iostream>
#include <iomanip>
#include <cassert>
#include <numcxx/numcxx.hxx>

int main()
{
    const int n=3;

    /// Declare four arrays with double elements.
    /// Like in the case with std::vector, the contents of a class
    /// object is placed on the stack, but the large data array is 
    /// put on the heap.

    numcxx::TArray1<double> A(n);
    numcxx::TArray1<double> B(n);
    numcxx::TArray1<double> C(n);
    //numcxx::TArray1<double> D(n);

    /// Declare and initialize  a matrix
    numcxx::TMatrix<double> Mat{{1,2,3},{4,5,6},{7,8,9}};

    // We can initialize the vectors with constant values.
    A=0.0;
    B=1.0;
    C=2.0;
    A=B+3.1;
    C=3.1+B;
    assert(numcxx::norm2(A-C)==0);
    A(1)=19;

    auto D = numcxx::arrayexpr(A*3);
//    std::cout << D;
    //D =A*3.0;
    std::cout << typeid(decltype(D)).name() << std::endl;
    std::cout << D[1] << std::endl;
    A=D*2.0;
//    A=3*(B+C/3.0)-(2.0*D);
    std::cout << std::endl<< std::endl;
    // Two ways to perform matrix-vector multiplication
    // Via expression template
    B=Mat*A;
    // Via apply method (calling dgemm of lapack)
    Mat.apply(A,C);

    
    std::cout<< B << std::endl;
    std::cout<< C << std::endl;

    auto diffnorm=numcxx::norm2(B-C);
    std::cout << "|B-C|="<<std::setprecision(10) <<  diffnorm << std::endl;
}
