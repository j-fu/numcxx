#include <cstdio>
#include <vector>

// initialize vector x with some data
void initialize(std::vector<double>& x)
{
    const int n=x.size();
    for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

// calculate the sum of the elements of x
double sum_elements(std::vector<double>& x)
{
    double sum=0;
    for (int i=0;i<x.size();i++)sum+=x[i];
    return sum;
}

int main()
{
    const int n=1.0e7;
    std::vector<double> x(n); // Construct vector with n elements
                              // Object "lives" on stack, data on heap
    initialize(x);
    double s=sum_elements(x);
    printf("sum=%e\n",s);
    // Object destructor automatically called at end of lifetime
    // So data array is freed automatically
}
