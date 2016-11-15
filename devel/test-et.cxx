#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <iostream>
#include "numcxx/numcxx.hxx"



int main()
{
    const int N=3;
    numcxx::DArray1 A(N);
    //numcxx::DArray1 B(N);
    numcxx::DArray1 C(N);
    numcxx::DArray1 D(N);

    numcxx::TMatrix<double> Mat{{1,2,3},{4,5,6},{7,8,9}};

    auto x=Mat;
    Mat(1,1)=100;
    std::cout <<Mat<< std::endl;
    std::cout <<x<< std::endl;
    A=0.0;
    auto B=A;
    B=1.0;
    std::cout<< B;
    C=2.0;
    D=3.25;
    A=B+3.1;
    std::cout<< B;
    std::cout<< A;
    A=3.1+B;
    std::cout<< A;
    D=A*3.0;
    A=3*(B+C/3.0)-4.0*D/2.0;

    std::cout<< A;
    
    B=3.0*(Mat*A);
    Mat.apply(A,C);
    std::cout<< B << std::endl;
    
   

    std::cout<< C << std::endl;

    // std::vector<double> E{33,44,55};
    // B=A+E;
    // std::cout<< B << std::endl;
    // B=Mat*E;
    // std::cout<< B << std::endl;
    // B=E;
    // std::cout<< B << std::endl;
}