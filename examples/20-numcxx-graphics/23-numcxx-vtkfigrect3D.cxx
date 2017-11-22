#include <numcxx/numcxx.hxx>

#include <cstdio>
#include <cmath>
#include "vtkfig/vtkfigFrame.h"
#include "vtkfig/vtkfigDataSet.h"
#include "vtkfig/vtkfigScalarView.h"


int main(void)
{
  const double x0=0.0;
  const double x1=1.0;
  double h=0.05;
  double t0=0.0;
  double t1=10.0;
  double dt=0.01;

  const int N=1+ceil((x1-x0)/h);
  h=(x1-x0)/(double)(N-1);

  numcxx::DArray1 X(N);
  numcxx::DArray1 U(N*N*N);

  X(0)=x0;
  for (int i=1;i<N;i++)
  {
    X(i)=X(i-1)+h;
  }
  
  

  auto griddata=vtkfig::DataSet::New();
  griddata->SetRectilinearGrid(X,X,X);

  
  auto frame=vtkfig::Frame::New();
  auto scview=vtkfig::ScalarView::New();
  scview->SetData(griddata,"U");
  frame->AddFigure(scview);
  
  double t=t0;


  while (t<t1)
  {
    int ijk=0;
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
        for (int k=0; k<N; k++,ijk++)
        {
          U(ijk)=sin(10.0*(X(i)-t))*cos(10.0*(X(j)-2*t))*cos(10.0*X(k)-0.5*t);
        }
    char titlebuf[20];
    snprintf(titlebuf,20,"t=%g",t);
    griddata->SetPointScalar(U ,"U");
    frame->Show();
    t+=dt;
  }
  
}
