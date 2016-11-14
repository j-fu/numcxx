#include <cstdio>
#include <vector>
#include <memory>

// initialize vector x with some data
void initialize(std::vector<double> &x)
{
    const int n=x.size();
    for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

// calculate the sum of the elements of x
double sum_elements(std::vector<double> & x)
{
    double sum=0;
    for (int i=0;i<x.size();i++)sum+=x[i];
    return sum;
}

int main()
{
    const int n=1.0e7;
    // call constructor and wrap pointer into smart pointer
    auto x=std::make_shared<std::vector<double>>(n);
    initialize(*x);  // dereference pointer to obtain reference
    double s=sum_elements(*x); // dereference pointer to obtain reference
    printf("sum=%e\n",s);
    // smartpointer calls desctrutor if reference count
    // reaches zero
}
