#include <cstdio>
#include <cstdlib>
#include <new>

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
    printf("xxx\n"); 
    const int n=1.0e7;
    double* x;
    try 
    {
        x=new double[n]; // allocate memory for vector on heap
    }
    catch (std::bad_alloc)
    {
        printf("error allocating x\n"); 
        exit(EXIT_FAILURE);
    }
    initialize(x,n);
    double s=sum_elements(x,n);
    printf("sum=%e\n",s);
    delete[] x; // Never forget to delete allocated memory
}
