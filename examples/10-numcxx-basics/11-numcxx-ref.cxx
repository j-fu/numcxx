///
/// \example numcxx-ref.cxx
///

#include <cstdio>
#include <numcxx/numcxx.hxx>

// initialize vector x with some data
void initialize(numcxx::DArray1 &X)
{
    const int n=X.size();
    for (int i=0;i<n;i++) X[i]= 1.0/(double)(1+n-i);
}

// calculate the sum of the elements of x
double sum_elements(numcxx::DArray1 & X)
{
    double sum=0;
    for (int i=0;i<X.size();i++)sum+=X[i];
    return sum;
}

int main()
{
    const int n=1.0e7;
    // call constructor 
    numcxx::TArray1<double> X(n);
    initialize(X);  
    double s=sum_elements(X);
    printf("sum=%e\n",s);
}
