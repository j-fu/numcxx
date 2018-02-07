///
/// \example 42-convtest-fem-sin.cxx
///
/// Stationary convergence test for fem method 
/// 
#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include <numcxx/numcxx.hxx>
#include <numcxx/simplegrid.hxx>
#include <numcxx/fem2d.hxx>
#include <netlib/netlib.hxx>

#ifdef VTKFIG
#include "numcxx/vtkfig-simplegrid.hxx"

#include "vtkfigFrame.h"
#include "vtkfigDataSet.h"
#include "vtkfigGridView.h"
#include "vtkfigScalarView.h"
#include "vtkfigXYPlot.h"
#endif

int main(void)
{
  double maxref=9;
  numcxx::Geometry Geometry;
  Geometry.set_points({
        {0,0},
        {1,0},
        {1,1},
        {0,1},
        {0.7,0.7},
    });

  Geometry.set_bfaces({
        {0,1},
        {1,2},
        {2,3},
        {3,0}
        });
  
  
  Geometry.set_bfaceregions({1,2,3,4});
  
  Geometry.set_regionpoints({
      {0.5,0.5}
    });
  Geometry.set_regionnumbers({1});


  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,800);
  frame->SetLayout(2,2);

  auto gridview=vtkfig::GridView::New();
  frame->AddFigure(gridview,0);

  auto solview=vtkfig::ScalarView::New();
  frame->AddFigure(solview,1);
  
  auto errplot=vtkfig::XYPlot::New();
  frame->AddFigure(errplot,2);
  errplot->SetXTitle("log10(h)");
  errplot->SetYTitle("log10(error)");

  auto timeplot=vtkfig::XYPlot::New();
  frame->AddFigure(timeplot,3);
  timeplot->SetXTitle("log10(N)");
  timeplot->SetYTitle("log10(time)");




  std::vector<double> H;
  std::vector<double> Cond;
  std::vector<double> N;
  std::vector<double> TAsm;
  std::vector<double> TGen;
  std::vector<double> TJac;
  std::vector<double> TCG;


  std::vector<double> L2Error;
  std::vector<double> H1Error;
  
  
  for (int iref=0;iref<maxref;iref++)
  {
    double t0;
    double h_intended=0.5*pow(2.0,-iref);
    double vol= h_intended*h_intended*0.25;
    
    Geometry.set_regionvolumes({vol});
    t0=numcxx::cpu_clock();
    numcxx::SimpleGrid grid(Geometry,"zpaAqDQ");
    double t_gen=numcxx::cpu_clock()-t0;

    double hmin,hmax;
    grid.calc_hminmax(hmin,hmax);
    printf("h_intended: %8.3g h_max: %8.3g\n",h_intended,hmax);
    
    
    numcxx::DArray1 bcfac(6);
    numcxx::DArray1 bcval(6);
    bcfac=0;
    bcval=0;
    
    bcfac(1)=fem2d::Dirichlet;
    bcfac(2)=fem2d::Dirichlet;
    bcfac(3)=fem2d::Dirichlet;
    bcfac(4)=fem2d::Dirichlet;
    bcval(1)=0.0;
    bcval(2)=0.0;
    bcval(3)=0.0;
    bcval(4)=0.0;
    
    auto nnodes=grid.npoints();
    auto nodes=grid.get_points();
    
    numcxx::DArray1 source(nnodes);
    numcxx::DArray1 exact(nnodes);
    numcxx::DArray1 kappa(nnodes);
    kappa=1;
    
    for (int i=0; i<nnodes;i++)
    {
      double x=nodes(i,0);
      double y=nodes(i,1);
      exact[i]=sin(M_PI*x)*sin(M_PI*y);
      source[i]=2.0*M_PI*M_PI*exact[i];
    }
    
    numcxx::DSparseMatrix SGlobal(nnodes,nnodes);
    numcxx::DArray1 Rhs(nnodes);
    numcxx::DArray1 Sol(nnodes);
    numcxx::DArray1 Res(nnodes);
    numcxx::DArray1 Upd(nnodes);
    
    t0=numcxx::cpu_clock();
    fem2d::assemble_heat_problem(grid,bcfac,bcval,source,kappa,SGlobal, Rhs);
    double t_asm=numcxx::cpu_clock()-t0;

    numcxx::DPreconJacobi Precon(SGlobal);

    t0=numcxx::cpu_clock();
    Precon.update(SGlobal);

    t0=numcxx::cpu_clock();
    Sol=0.0;
    unsigned long max_iter=10000000;
    double eps=1.0e-10;
    double res0=numcxx::norm2(SGlobal*Sol-Rhs);
    unsigned long nsteps=0;
    for (;nsteps<max_iter;nsteps++)
    {
      Res=SGlobal*Sol-Rhs;
      Precon.solve(Upd,Res);
      Sol=Sol-Upd;
      double norm=numcxx::norm2(Upd);
      if (norm<1.0e-10) break;
    }
    double res=numcxx::norm2(SGlobal*Sol-Rhs);
    double rho=pow(res/res0,1.0/nsteps);
    double cond=(1.0+rho)/(1.0-rho);
    printf("Jacobi contr: %g, per step: %g, kappa: %g\n", res/res0, rho,cond);
    double t_jac=numcxx::cpu_clock()-t0;



    t0=numcxx::cpu_clock();
    Sol=0.0;
    res0=numcxx::norm2(SGlobal*Sol-Rhs);
    int max_iter_cg=100000;
    netlib::CG(SGlobal,Sol, Rhs,Precon,max_iter_cg,eps);
    res=numcxx::norm2(SGlobal*Sol-Rhs);
    rho=pow(res/res0,1.0/max_iter_cg);
    printf("CG contr: %g, per step: %g\n", res/res0, rho);
    double t_cg=numcxx::cpu_clock()-t0;


    



    double l2error=fem2d::l2norm(grid,Sol-exact);
    double h1error=fem2d::h1norm(grid,Sol-exact);

    printf("h: %8.3e l2: %8.3e h1: %8.3e\n",hmax,l2error,h1error);
    H.push_back(log10(hmax));
    H1Error.push_back(log10(h1error));
    L2Error.push_back(log10(l2error));
    Cond.push_back(log10(cond));

    printf("time/ms gen: %8.3f asm: %8.3f jac: %8.3f cg: %8.3f\n",
           t_gen*1000.0, t_asm*1000.0, t_jac*1000.0,t_cg*1000.0);
    
    N.push_back(log10(grid.npoints()));;
    TAsm.push_back(log10(t_asm));
    TJac.push_back(log10(t_jac));
    TCG.push_back(log10(t_cg));
    TGen.push_back(log10(t_gen));


#ifdef VTKFIG
    auto griddata=numcxx::vtkfigDataSet(grid);
    griddata->SetPointScalar(Sol ,"Sol");

    gridview->SetData(griddata);
    solview->SetData(griddata,"Sol");

    frame->SetFrameTitle("iref="+std::to_string(iref));
    if (iref>0)
    {
      double x0=-3; 
      double x1=0;
      auto X=numcxx::DArray1{x0,x1};
      auto OH1=numcxx::DArray1{x0,x1};
      auto OH2=numcxx::DArray1{2*x0,2*x1};
      auto OHM2=numcxx::DArray1{-2*x0-1,-2*x1-1};
      
      errplot->Clear();


      errplot->SetLegendPosition(0.3,0.4);
      errplot->SetXRange(x0,x1);
      errplot->SetYRange(-5,5);


      errplot->SetPlotLineType("-");
      errplot->SetPlotMarkerType("o");
      errplot->SetMarkerSize(2);


      errplot->SetPlotColor(0.2,0.5,0.2);
      errplot->SetPlotLegend("O(h**2)");
      errplot->AddPlot(X,OH2);
      errplot->SetPlotColor(0.5,0.2,0.2);
      errplot->SetPlotLegend("O(h)");
      errplot->AddPlot(X,OH1);

      errplot->SetPlotColor(0.0,0.5,0);
      errplot->SetPlotLegend("L2");
      errplot->AddPlot(H,L2Error);

      errplot->SetPlotLegend("H1");
      errplot->SetPlotColor(0.5,0,0);
      errplot->AddPlot(H,H1Error);

      errplot->SetPlotLegend("Cond");
      errplot->SetPlotColor(0,0,0.5);
      errplot->AddPlot(H,Cond);
      errplot->SetPlotLegend("O(1/h**2)");
      errplot->AddPlot(X,OHM2);

      
      double n0=0;
      double n1=7;
      double n12=4;
      double n11=6;
      double t0=-5;
      double t1=2;

      auto XN=numcxx::DArray1{n0,n1};
      auto XN1=numcxx::DArray1{n0,n11};
      auto XN2=numcxx::DArray1{n0,n12};
      auto ON1=numcxx::DArray1{n0-6.0,n1-6.0};
      auto ON2=numcxx::DArray1{2.0*n0-6, 2*n12-6};
      auto ON32=numcxx::DArray1{1.5*n0-8.0, 1.5*n11-8.0};


      timeplot->SetLegendPosition(0.3,0.6);
      timeplot->SetLegendSize(0.15,0.3);
      timeplot->SetXRange(n0,n1);
      timeplot->SetYRange(t0,t1);


      timeplot->Clear();
      timeplot->SetPlotLineType("-");
      timeplot->SetPlotMarkerType("o");
      timeplot->SetMarkerSize(2);


      // timeplot->SetPlotLegend("Gen");
      // timeplot->SetPlotColor(0.5,0,0);
      // timeplot->AddPlot(N,TGen);

      // timeplot->SetPlotLegend("Asm");
      // timeplot->SetPlotColor(0.0,0.5,0);
      // timeplot->AddPlot(N,TAsm);

      timeplot->SetPlotLegend("Jacobi");
      timeplot->SetPlotColor(0.0,0.0,0.5);
      timeplot->AddPlot(N,TJac);

      timeplot->SetPlotLegend("CG");
      timeplot->SetPlotColor(0.0,0.5,0.5);
      timeplot->AddPlot(N,TCG);

      timeplot->SetPlotLegend("O(N**2))");
      timeplot->SetPlotColor(0.0,0.0,1.0);
      timeplot->AddPlot(XN2,ON2);

      timeplot->SetPlotLegend("O(N**(3/2))");
      timeplot->SetPlotColor(1.0,0.0,0.0);
      timeplot->AddPlot(XN1,ON32);

      // timeplot->SetPlotLegend("O(N)");
      // timeplot->SetPlotColor(1.0,0,0);
      // timeplot->AddPlot(XN,ON1);;


    }
    frame->Interact();
#endif
  }

}


