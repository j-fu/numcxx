#include <iostream>
#include <iomanip>
#include <limits>
#include <ctgmath>


template <typename T> void basel(long n)
{
    const T PI=3.1415926535897932384626433L;

    // calculate sum in forward direction
    T sum_fwd=0.0;
    for (int i=1; i<=n;i++)
    {   T x=i;
        sum_fwd += 1.0/(x*x);

    }

    // calculate sum in reverse direction: 
    // first adding up small members leads to increased
    // precision
    T sum_rev=0.0;
    for (int i=n; i>=1;i--)
    {   T x=i;
        sum_rev += 1.0/(x*x);
    }



    T sum_exact=PI*PI/6.0;
    std::cout 
        << std::scientific
        << std::setprecision(20)
        << std::setw(20)
        << n 
        << " " 
        << sum_fwd 
        << " " 
        << std::fabs(sum_exact-sum_fwd)
        << " " 
        << sum_rev
        << " " 
        << std::fabs(sum_exact-sum_rev)
        << std::endl;
}

int main()
{
    std::cout << "Float version:" << std::endl;
    std::cout << std::setw(20)
              << "n"
              << " " 
              << std::setw(24)
              << "forward sum" 
              << " " 
              << std::setw(24)
              << "forward sum error"
              << " " 
              << std::setw(24)
              << "reverse sum"
              << " " 
              << std::setw(24)
              << "reverse sum error"
              << std::endl;
    
    for (int i=0,n=10;i<8;i++,n*=10) basel<float>(n);
    std::cout << "Double version:" << std::endl;
    for (int i=0,n=10;i<8;i++,n*=10) basel<double>(n);
}
