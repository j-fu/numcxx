#include <iostream>
#include <iomanip>
#include <limits>


template <typename T> 
T print_macheps(void)
{
    T eps=1.0;
    T one=1.0;
    T epsnew=1.0;
    T result=0.0;

    do 
    {
        eps=epsnew;
        epsnew=eps/2.0;
        result=one+epsnew;
    } while (result>one);
    
    // we need to be careful here: 
    // the value of epsnew at which one+epsnew = one
    // is half of the value we searched for...

    std::cout << "   Calculated: "
              << std::setprecision(30)
              << eps 
              << std::endl
              << "From <limits>: " 
              << std::numeric_limits<T>::epsilon() 
              << std::endl
              << std::endl;

    return eps;
}


int main(int argc, const char *argv[])
{
    print_macheps<float>();
    print_macheps<double>();
    print_macheps<long double>();
    return 0;
}
