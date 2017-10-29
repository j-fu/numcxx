#include <iostream>
#include <numcxx/numcxx.hxx>

/// This is a test for the expression templates implemented
/// in numcxx. 
int main()
{
    const int n=3;
    numcxx::DArray1 A(n);
    numcxx::DArray1 B(n);
    numcxx::DArray1 C(n);
    numcxx::DArray1 D(n);

    numcxx::TMatrix<double> Mat{{1,2,3},{4,5,6},{7,8,9}};
    A=0.0;
    B=1.0;
    C=2.0;
    D=3.25;
    A=B+3.1;
    A=3.1+B;
    D=A*3.0;
    A=3*(B+C/3.0)-4.0*D/2.0;
    
    // Two ways to perform matrix-vector multiplication
    // Via expression template
    B=Mat*A;
    // Via apply method (calling dgemm of lapack)
    Mat.apply(A,C);

    std::cout<< B << std::endl;
    std::cout<< C << std::endl;
}
