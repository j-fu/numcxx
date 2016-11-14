#include <cstdio>

// initialize vector x with some data
void initialize(double *x, int n)
{
    for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

// calculate the sum of the elements of x
double sum_elements(double *x, int n)
{
    double sum=0;
    for (int i=0;i<n;i++) sum+=x[i];
    return sum;
}


int main()
{
    const int n=1.0e7; // double value automatically converted to int
    double x[n];       // declare array on stack
    initialize(x,n);
    double s=sum_elements(x,n);
    printf("sum=%e\n",s);
    // no need to free variable defined on stack
}
