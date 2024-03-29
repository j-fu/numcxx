#include <numcxx/numcxx.hxx>

#include <cstdio>
#include <cmath>
#include "vtkfig/vtkfigFrame.h"
#include "vtkfig/vtkfigXYPlot.h"


int main(void)
{
  const double x0=0.0;
  const double x1=1.0;
  double h=0.01;
  double t0=0.0;
  double t1=10.0;
  double dt=0.01;

  const int N=1+ceil((x1-x0)/h);
  h=(x1-x0)/(double)(N-1);

  numcxx::DArray1 X(N);
  numcxx::DArray1 U(N);

  X(0)=x0;
  for (int i=1;i<N;i++)
  {
    X(i)=X(i-1)+h;
  }
  
  
  
  auto frame=vtkfig::Frame::New();
  auto xyplot=vtkfig::XYPlot::New();
  xyplot->SetPlotColor(1.0,0.0,0.0);
  xyplot->SetPlotLineType("-");
  xyplot->SetXAxisLabelFormat("%3.1f");
  xyplot->SetYAxisLabelFormat("%3.1f");
  frame->AddFigure(xyplot);
  
  double t=t0;


  while (t<t1)
  {
    for (int i=0; i<N; i++)
    {
      U(i)=sin(20.0*(X(i)-t));
    }
    
    char titlebuf[20];
    snprintf(titlebuf,20,"t=%g",t);
    xyplot->Clear();
    xyplot->SetTitle(titlebuf);
    xyplot->AddPlot(X,U);
    frame->Show();
    t+=dt;
  }
  
}


