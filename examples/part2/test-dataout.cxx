#include <cstdio>
#include <cmath>
#include <numcxx/numcxx.hxx>



int main()
{
    const int n=21;
    numcxx::TArray1<double> X(n),Y(n);


    double L=10.0;
    double h=L/(double) (n-1);
    
    double x=0.0;    
    for (int i=0;i<n;i++,x+=h)
        X(i)=x;

    for (int i=0;i<n;i++)
        Y(i)=sin(X(i));

    
    FILE* out=fopen("test.dat","w");
    for (int i=0;i<n;i++)
    {
        fprintf(out,"%e %e\n",X(i),Y(i));
    }
    fclose(out);

}
