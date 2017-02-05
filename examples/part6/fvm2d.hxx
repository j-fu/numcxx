#ifndef FEM2D_H
#define FEM2D_H

#include <numcxx/simplegrid.hxx>
#include <numcxx/tsparsematrix.hxx>
#include <cmath>
#include <iostream>

namespace fvm2d
{
  // Assemble stiffness matrix for mixed Dirichlet/Homogeneus Neumann
  // heat conduction
  const double Dirichlet=1.0e30;

  inline void compute_local_formfactors(
    const int icell,
    const numcxx::DArray2 & points,
    const numcxx::IArray2 & cells,
    numcxx::DArray1 & epar,
    numcxx::DArray1 & npar,
    double & vol)
  {
    int i0=cells(icell,0);
    int i1=cells(icell,1);
    int i2=cells(icell,2);
    
    
    // Fill matrix of edge vectors
    double V00= points(i1,0)- points(i0,0);
    double V10= points(i1,1)- points(i0,1);
    
    double V01= points(i2,0)- points(i0,0);
    double V11= points(i2,1)- points(i0,1);
    
    double V02= points(i2,0)- points(i1,0);
    double V12= points(i2,1)- points(i1,1);
    
    
    
    // Compute determinant 
    double det=V00*V11 - V01*V10;
    vol=0.5*det;
    
    double ivol = 1.0/vol;
    
    // squares of edge lengths
    double dd0=V02*V02+V12*V12; // l21
    double dd1=V01*V01+V11*V11; // l20
    double dd2=V00*V00+V10*V10; // l10
    
    
    // contributions to \sigma_kl/h_kl
    epar(0)= (dd1+dd2-dd0)*0.125*ivol;
    epar(1)= (dd2+dd0-dd1)*0.125*ivol;
    epar(2)= (dd0+dd1-dd2)*0.125*ivol;
    
    
    // contributions to \omega_k
    npar(0)= (epar(2)*dd2+epar(1)*dd1)*0.25;
    npar(1)= (epar(0)*dd0+epar(2)*dd2)*0.25;
    npar(2)= (epar(1)*dd1+epar(0)*dd0)*0.25;
  }
  


  inline void
  assemble_heat_problem_with_source(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSource,
    std::shared_ptr<numcxx::DArray1> pKappa,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    )
    {
    auto & grid=*pGrid;
    auto & bcfac=*pBCfac;
    auto & bcval=*pBCval;
    auto & Rhs=*pRhs;
    auto & S=*pS;
    auto & source=*pSource;
    auto & kappa=*pKappa;

    auto ndim=grid.spacedim();

    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();


    auto pEpar=numcxx::DArray1::create(ndim+1);
    auto pNpar=numcxx::DArray1::create(ndim+1);
    
    auto &epar=*pEpar; // s_i/h_i
    auto &npar=*pNpar; // omega_i
    double vol;

    Rhs=0.0;
    S(0,0)=0;
    S.clear();

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
           
      int k0=cells(icell,0);
      int k1=cells(icell,1);
      int k2=cells(icell,2);
      compute_local_formfactors(icell, points,cells,epar,npar, vol);

      Rhs(k0)+=source(k0)*npar(0);
      Rhs(k1)+=source(k1)*npar(1);
      Rhs(k2)+=source(k2)*npar(2);

      // Assemble fluxes with edge averaged diffusion
      // coefficients into global matrix
      S(k0,k0)+=epar(2)*0.5*(kappa(k0)+kappa(k1));
      S(k0,k1)-=epar(2)*0.5*(kappa(k0)+kappa(k1));
      S(k1,k0)-=epar(2)*0.5*(kappa(k0)+kappa(k1));
      S(k1,k1)+=epar(2)*0.5*(kappa(k0)+kappa(k1));

      S(k0,k0)+=epar(1)*0.5*(kappa(k0)+kappa(k2));
      S(k0,k2)-=epar(1)*0.5*(kappa(k0)+kappa(k2));
      S(k2,k0)-=epar(1)*0.5*(kappa(k0)+kappa(k2));
      S(k2,k2)+=epar(1)*0.5*(kappa(k0)+kappa(k2));

      S(k1,k1)+=epar(0)*0.5*(kappa(k1)+kappa(k2));
      S(k1,k2)-=epar(0)*0.5*(kappa(k1)+kappa(k2));
      S(k2,k1)-=epar(0)*0.5*(kappa(k1)+kappa(k2));
      S(k2,k2)+=epar(0)*0.5*(kappa(k1)+kappa(k2));
    }    


    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
      double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
      double h=sqrt(dx*dx+dy*dy);
      double factor=0.5*h*bcfac(ireg);
      if (bcfac(ireg)>=Dirichlet) factor=bcfac(ireg);

      int i;
      i=bfaces(ibface,0);
      S(i,i)+=factor;
      Rhs(i)+=factor*bcval(ireg);

      i=bfaces(ibface,1);
      S(i,i)+=factor;
      Rhs(i)+=factor*bcval(ireg);
    }
    S.flush();
  }




/////////////////////////////////////////////////

  inline double kappa(double u)
  {
    return 0.1+100*u*u*u*u;
  }
  inline double dkappa(double u)
  {
    return 400.0*u*u*u; 
  }


  void assemble_and_apply_nonlinear_heat(
    std::shared_ptr<numcxx::SimpleGrid >pGrid,// Discretization grid
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSol,
    std::shared_ptr<numcxx::DSparseMatrix> pS,
    std::shared_ptr<numcxx::DArray1> pRhs
    )
  {
    auto & grid=*pGrid;
    auto & bcfac=*pBCfac;
    auto & bcval=*pBCval;
    auto & Rhs=*pRhs;
    auto & S=*pS;
    auto & Sol=*pSol;

    auto ndim=grid.spacedim();

    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();
    Rhs=0.0;
    S(0,0)=0;
    S.clear();



    auto pEpar=numcxx::DArray1::create(ndim+1);
    auto pNpar=numcxx::DArray1::create(ndim+1);
    
    auto &epar=*pEpar;
    auto &npar=*pNpar;
    
    double vol;

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
      
      compute_local_formfactors(icell, points,cells,epar,npar, vol);

      int k0=cells(icell,0);
      int k1=cells(icell,1);
      int k2=cells(icell,2);

      double kappa0=kappa(Sol(k0));
      double kappa1=kappa(Sol(k1));
      double kappa2=kappa(Sol(k2));

      double dkappa0=dkappa(Sol(k0));
      double dkappa1=dkappa(Sol(k1));
      double dkappa2=dkappa(Sol(k2));
      

      double flux0=epar(0)*0.5*(kappa1+kappa2)*(Sol(k1)-Sol(k2));
      Rhs(k1)+=flux0;
      Rhs(k2)-=flux0;
      S(k1,k1)+= epar(0)*( 0.5*(kappa1+kappa2)+0.5*dkappa1*(Sol(k1)-Sol(k2)));
      S(k1,k2)+= epar(0)*(-0.5*(kappa1+kappa2)+0.5*dkappa2*(Sol(k1)-Sol(k2)));
      S(k2,k1)+= epar(0)*(-0.5*(kappa1+kappa2)-0.5*dkappa1*(Sol(k1)-Sol(k2)));
      S(k2,k2)+= epar(0)*( 0.5*(kappa1+kappa2)-0.5*dkappa2*(Sol(k1)-Sol(k2)));

      double flux1=epar(1)*0.5*(kappa2+kappa0)*(Sol(k2)-Sol(k0));
      Rhs(k2)+=flux1;
      Rhs(k0)-=flux1;
      S(k2,k2)+= epar(1)*( 0.5*(kappa0+kappa2)+0.5*dkappa2*(Sol(k2)-Sol(k0)));
      S(k2,k0)+= epar(1)*(-0.5*(kappa0+kappa2)+0.5*dkappa0*(Sol(k2)-Sol(k0)));
      S(k0,k2)+= epar(1)*(-0.5*(kappa0+kappa2)-0.5*dkappa2*(Sol(k2)-Sol(k0)));
      S(k0,k0)+= epar(1)*( 0.5*(kappa0+kappa2)-0.5*dkappa0*(Sol(k2)-Sol(k0)));


      double flux2=epar(2)*0.5*(kappa0+kappa1)*(Sol(k0)-Sol(k1));
      Rhs(k0)+=flux2;
      Rhs(k1)-=flux2;

      S(k0,k0)+= epar(2)*( 0.5*(kappa0+kappa1)+0.5*dkappa0*(Sol(k0)-Sol(k1)));
      S(k0,k1)+= epar(2)*(-0.5*(kappa0+kappa1)+0.5*dkappa1*(Sol(k0)-Sol(k1)));
      S(k1,k0)+= epar(2)*(-0.5*(kappa0+kappa1)-0.5*dkappa0*(Sol(k0)-Sol(k1)));
      S(k1,k1)+= epar(2)*( 0.5*(kappa0+kappa1)-0.5*dkappa1*(Sol(k0)-Sol(k1)));
    }    


    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
      double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
      double h=sqrt(dx*dx+dy*dy);
      double factor=0.5*h*bcfac(ireg);
      if (bcfac(ireg)>=Dirichlet) factor=bcfac(ireg);

      int i;
      i=bfaces(ibface,0);
      S(i,i)+=factor;
      Rhs(i)+=factor*(Sol(i)-bcval(ireg));

      i=bfaces(ibface,1);
      S(i,i)+=factor;
      Rhs(i)+=factor*(Sol(i)-bcval(ireg));
    }
    
    S.flush();
  }

  void initialize_bc(
    numcxx::SimpleGrid &grid,// Discretization grid
    numcxx::DArray1& g,
    numcxx::DArray1& Sol
    )
  {

    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int i;
      int ireg=bfaceregions(ibface);


      i=bfaces(ibface,0);
      Sol(i)=g(ireg);

      i=bfaces(ibface,1);
      Sol(i)=g(ireg);
    }

  }

////////////////////////////////////////////////////////






  void initialize_bc(
    std::shared_ptr<numcxx::SimpleGrid >grid,// Discretization grid
    std::shared_ptr<numcxx::DArray1> g,
    std::shared_ptr<numcxx::DArray1> Sol
    )
  {
    initialize_bc(*grid,*g, *Sol);
  }

}
#endif
