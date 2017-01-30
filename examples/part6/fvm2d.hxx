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
  const double DirichletPenalty=1.0e30;

 
  // Assemble stiffness matrix for mixed Dirichlet/Homogeneus Neumann
  // heat conduction
  inline std::shared_ptr<numcxx::TSparseMatrix<double>> 
  assemble_general_heat_matrix(
    numcxx::SimpleGrid &grid,  // Discretization grid
    numcxx::DArray1& kappa, // heat conduction coefficient (per node)
    numcxx::DArray1& alpha // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    )
  {
    auto ndim=grid.spacedim();

    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();


    // Global stiffness matrix
    auto pSGlobal=numcxx::DSparseMatrix::create(npoints, npoints);
    auto &SGlobal=*pSGlobal;

    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto &V=*pV;

    auto pEpar=numcxx::DArray1::create(ndim+1);
    auto pNpar=numcxx::DArray1::create(ndim+1);
    auto pDD=numcxx::DArray1::create(ndim+1);
    
    auto &epar=*pEpar;
    auto &npar=*pNpar;
    auto &dd=*pDD;

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
           
      
      // Fill matrix of edge vectors
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);

      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);

      V(0,2)= points(cells(icell,2),0)- points(cells(icell,1),0);
      V(1,2)= points(cells(icell,2),1)- points(cells(icell,1),1);

      

      // Compute determinant 
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double ivol = 2.0/det;

      // squares of edge lengths
      dd(0)=V(0,2)*V(0,2)+V(1,2)*V(1,2); // l21
      dd(1)=V(0,1)*V(0,1)+V(1,1)*V(1,1); // l20
      dd(2)=V(0,0)*V(0,0)+V(1,0)*V(1,0); // l10

      
      // contributions to \sigma_kl/h_kl
      epar(0)= (dd(1)+dd(2)-dd(0))*0.125*ivol;
      epar(1)= (dd(2)+dd(0)-dd(1))*0.125*ivol;
      epar(2)= (dd(0)+dd(1)-dd(2))*0.125*ivol;


      // contributions to \omega_k
      npar(0)= (epar(2)*dd(2)+epar(1)*dd(1))*0.25;
      npar(1)= (epar(0)*dd(0)+epar(2)*dd(2))*0.25;
      npar(2)= (epar(1)*dd(1)+epar(0)*dd(0))*0.25;

      int k0=cells(icell,0);
      int k1=cells(icell,1);
      int k2=cells(icell,2);

      // Assemble fluxes with edge averaged diffusion
      // coefficients into global matrix
      SGlobal(k0,k0)+=epar(2)*0.5*(kappa(k0)+kappa(k1));
      SGlobal(k0,k1)-=epar(2)*0.5*(kappa(k0)+kappa(k1));
      SGlobal(k1,k0)-=epar(2)*0.5*(kappa(k0)+kappa(k1));
      SGlobal(k1,k1)+=epar(2)*0.5*(kappa(k0)+kappa(k1));

      SGlobal(k0,k0)+=epar(1)*0.5*(kappa(k0)+kappa(k2));
      SGlobal(k0,k2)-=epar(1)*0.5*(kappa(k0)+kappa(k2));
      SGlobal(k2,k0)-=epar(1)*0.5*(kappa(k0)+kappa(k2));
      SGlobal(k2,k2)+=epar(1)*0.5*(kappa(k0)+kappa(k2));

      SGlobal(k1,k1)+=epar(0)*0.5*(kappa(k1)+kappa(k2));
      SGlobal(k1,k2)-=epar(0)*0.5*(kappa(k1)+kappa(k2));
      SGlobal(k2,k1)-=epar(0)*0.5*(kappa(k1)+kappa(k2));
      SGlobal(k2,k2)+=epar(0)*0.5*(kappa(k1)+kappa(k2));
    }    


    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);

      // Check if it is "Dirichlet"
      if (alpha(ireg)>=DirichletPenalty)
      {
        // Assemble penalty values
        int i;
        i=bfaces(ibface,0);
        SGlobal(i,i)+=alpha(ireg);

        i=bfaces(ibface,1);
        SGlobal(i,i)+=alpha(ireg);
      }
      else
      {
        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        // first oder quadrature, needs to be replaced by second order
        int i0=bfaces(ibface,0);
        int i1=bfaces(ibface,1);

        SGlobal(i0,i0)+=h*alpha(ireg)/2.0;
        SGlobal(i1,i1)+=h*alpha(ireg)/2.0;
      }
    }
    // Return global stiffness matrix
    return pSGlobal;
  }


  // Assemble zero right hand side for mixed Dirichlet/Homogeneus Neumann problem
  std::shared_ptr<numcxx::TArray1<double>>
  assemble_heat_rhs_zero(
    numcxx::SimpleGrid &g,       // Discretization grid
    numcxx::DArray1& alpha, // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    numcxx::DArray1 &bcval       // Dirichlet boundary values
    )
  {
    int npoints=g.npoints();
    auto points=g.get_points(); // Array of global nodes
    int nbfaces=g.nbfaces();
    auto bfaces=g.get_bfaces();
    auto bfaceregions=g.get_bfaceregions();

    auto pRhs=numcxx::DArray1::create(npoints);
    auto &Rhs=*pRhs;

    Rhs=0.0;

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int ireg=bfaceregions(ibface);
      if (alpha(ireg)>=DirichletPenalty)
      {
        int i;

        i=bfaces(ibface,0);
        Rhs(i)+=DirichletPenalty*bcval(ireg);

        i=bfaces(ibface,1);
        Rhs(i)+=DirichletPenalty*bcval(ireg);
      }
      else
      {
        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        // first oder quadrature, needs to be replaced by second order
        int i0=bfaces(ibface,0);
        int i1=bfaces(ibface,1);
        
        Rhs(i0)+=h*alpha(ireg)*bcval(ireg)/2.0;
        Rhs(i1)+=h*alpha(ireg)*bcval(ireg)/2.0;
      }
    }
    return pRhs;
  }
  


  inline double kappa(double u)
  {
   return 0.1+10*u*u;
  }
  inline double dkappa(double u)
  {
    return 20.0*u;
  }

  void assemble_and_apply_nonlinear_heat(
    numcxx::SimpleGrid &grid,// Discretization grid
    numcxx::DSparseMatrix& S,
    numcxx::DArray1& Sol,
    numcxx::DArray1& Rhs,
    numcxx::DArray1& alpha,
    numcxx::DArray1& g
    )
  {
    auto ndim=grid.spacedim();

    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();
    Rhs=0.0;
    S(0,0)=0;
    S.clear();


    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto &V=*pV;

    auto pEpar=numcxx::DArray1::create(ndim+1);
    auto pNpar=numcxx::DArray1::create(ndim+1);
    auto pDD=numcxx::DArray1::create(ndim+1);
    
    auto &epar=*pEpar;
    auto &npar=*pNpar;
    auto &dd=*pDD;

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
      
      // Fill matrix of edge vectors
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);

      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);

      V(0,2)= points(cells(icell,2),0)- points(cells(icell,1),0);
      V(1,2)= points(cells(icell,2),1)- points(cells(icell,1),1);


      // Compute determinant 
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double ivol = 2.0/det;

      dd(0)=V(0,2)*V(0,2)+V(1,2)*V(1,2); // l21
      dd(1)=V(0,1)*V(0,1)+V(1,1)*V(1,1); // l20
      dd(2)=V(0,0)*V(0,0)+V(1,0)*V(1,0); // l10

// 012

      epar(0)= (dd(1)+dd(2)-dd(0))*0.125*ivol;
      epar(1)= (dd(2)+dd(0)-dd(1))*0.125*ivol;
      epar(2)= (dd(0)+dd(1)-dd(2))*0.125*ivol;


      npar(0)= (epar(2)*dd(2)+epar(1)*dd(1))*0.25;
      npar(1)= (epar(0)*dd(0)+epar(2)*dd(2))*0.25;
      npar(2)= (epar(1)*dd(1)+epar(0)*dd(0))*0.25;


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
      double factor=0.5*h*alpha(ireg);
      if (alpha(ireg)>=DirichletPenalty) factor=alpha(ireg);

      int i;
      i=bfaces(ibface,0);
      S(i,i)+=factor;
      Rhs(i)+=factor*(Sol(i)-g(ireg));

      i=bfaces(ibface,1);
      S(i,i)+=factor;
      Rhs(i)+=factor*(Sol(i)-g(ireg));
    }
    
    S.flush();
  }

  void initialize_bc(
    numcxx::SimpleGrid &grid,// Discretization grid
    numcxx::DArray1& Sol,
    numcxx::DArray1& g
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


  inline 
  std::shared_ptr<numcxx::DSparseMatrix> 
  assemble_general_heat_matrix(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> kappa, // heat conduction coefficient (per node)
    std::shared_ptr<numcxx::DArray1> alpha // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    )
  {
    return assemble_general_heat_matrix(*grid,*kappa,*alpha);
  }


  // Assemble zero right hand side for mixed Dirichlet/Homogeneus Neumann problem
  std::shared_ptr<numcxx::DArray1>
  assemble_heat_rhs_zero(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> alpha,     // boundary heat transfer coefficient (large value marks Dirichlet)
    std::shared_ptr<numcxx::DArray1> g        // boundary ambient temperature
    )
  {
    return assemble_heat_rhs_zero(*grid,*alpha,*g);
  }



  void assemble_and_apply_nonlinear_heat(
    std::shared_ptr<numcxx::SimpleGrid >grid,// Discretization grid
    std::shared_ptr<numcxx::DSparseMatrix> S,
    std::shared_ptr<numcxx::DArray1> Sol,
    std::shared_ptr<numcxx::DArray1> Rhs,
    std::shared_ptr<numcxx::DArray1> alpha,
    std::shared_ptr<numcxx::DArray1> g
    )
  {
    assemble_and_apply_nonlinear_heat(*grid,*S,*Sol,*Rhs,*alpha,*g);
  }


  void initialize_bc(
    std::shared_ptr<numcxx::SimpleGrid >grid,// Discretization grid
    std::shared_ptr<numcxx::DArray1> Sol,
    std::shared_ptr<numcxx::DArray1> g
    )
  {
    initialize_bc(*grid,*Sol,*g);
  }

}
#endif
