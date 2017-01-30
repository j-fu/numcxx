#ifndef FEM2D_H
#define FEM2D_H

#include <numcxx/simplegrid.hxx>
#include <numcxx/tsparsematrix.hxx>
#include <cmath>
#include <iostream>

namespace fem2d
{
  // Assemble stiffness matrix for mixed Dirichlet/Homogeneus Neumann
  // heat conduction
  const double DirichletPenalty=1.0e30;

  inline std::shared_ptr<numcxx::TSparseMatrix<double>> 
  assemble_heat_matrix(
    numcxx::SimpleGrid &g,       // Discretization grid
    numcxx::IArray1 &bcdirichlet // Boundary regions with Dirichlet BC
    )
  {
    auto ndim=g.spacedim();

    auto points=g.get_points(); // Array of global nodes
    auto cells=g.get_cells();   // Local-global dof map

    int npoints=g.npoints();
    int ncells=g.ncells();


    // Global stiffness matrix
    auto pSGlobal=numcxx::DSparseMatrix::create(npoints, npoints);
    auto &SGlobal=*pSGlobal;

    // Local stiffness matrix
    auto pSLocal=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;

    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim, ndim);
    auto &V=*pV;

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
      
      // Fill matrix V
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);

      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);


      // Compute determinant
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double invvol = 2.0/det;

      
      // Compute entris of local stiffness matrix
      SLocal(0,0)= invvol * (  ( V(1,0)-V(1,1) )*( V(1,0)-V(1,1) )+( V(0,1)-V(0,0) )*( V(0,1)-V(0,0) ) );
      SLocal(0,1)= invvol * (  ( V(1,0)-V(1,1) )* V(1,1)          - ( V(0,1)-V(0,0) )*V(0,1) );
      SLocal(0,2)= invvol * ( -( V(1,0)-V(1,1) )* V(1,0)          + ( V(0,1)-V(0,0) )*V(0,0) );

      SLocal(1,1)=  invvol * (  V(1,1)*V(1,1) + V(0,1)*V(0,1) );
      SLocal(1,2)=  invvol * ( -V(1,1)*V(1,0) - V(0,1)*V(0,0) );

      SLocal(2,2)=  invvol * ( V(1,0)*V(1,0)+ V(0,0)*V(0,0) );

      SLocal(1,0)=SLocal(0,1);
      SLocal(2,0)=SLocal(0,2);
      SLocal(2,1)=SLocal(1,2);

      // Assemble into global stiffness matrix
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
          SGlobal(cells(icell,i),cells(icell,j))+=SLocal(i,j);
    }    


    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=g.nbfaces();
    auto bfaces=g.get_bfaces();
    auto bfaceregions=g.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);

      // Check if it is "Dirichlet"
      if (bcdirichlet(ireg))
      {
        // Assemble penalty values
        int i;
        i=bfaces(ibface,0);
        SGlobal(i,i)+=DirichletPenalty;

        i=bfaces(ibface,1);
        SGlobal(i,i)+=DirichletPenalty;
      }
    }
    // Return global stiffness matrix
    return pSGlobal;
  }

  // Assemble zero right hand side for mixed Dirichlet/Homogeneus Neumann problem
  std::shared_ptr<numcxx::TArray1<double>>
  assemble_heat_rhs(
    numcxx::SimpleGrid &g,       // Discretization grid
    numcxx::IArray1 &bcdirichlet,// Boundary regions with Dirichlet BC
    numcxx::DArray1 &bcval       // Dirichlet boundary values
    )
  {
    int npoints=g.npoints();
    int nbfaces=g.nbfaces();
    auto bfaces=g.get_bfaces();
    auto bfaceregions=g.get_bfaceregions();

    auto pRhs=numcxx::DArray1::create(npoints);
    auto &Rhs=*pRhs;

    Rhs=0.0;

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int ireg=bfaceregions(ibface);
      if (bcdirichlet(ireg))
      {
        int i;

        i=bfaces(ibface,0);
        Rhs(i)+=DirichletPenalty*bcval(ireg);

        i=bfaces(ibface,1);
        Rhs(i)+=DirichletPenalty*bcval(ireg);
      }
    }
    return pRhs;
  }
  

  /////////////////////////////////////////////////////////////////////////////////////

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

    // Local stiffness matrix
    auto pSLocal=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;

    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim, ndim);
    auto &V=*pV;

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
      
      // Fill matrix V
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);

      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);


      // Compute determinant
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double invvol = 2.0/det;

      
      // Compute entris of local stiffness matrix
      SLocal(0,0)= invvol * (  ( V(1,0)-V(1,1) )*( V(1,0)-V(1,1) )+( V(0,1)-V(0,0) )*( V(0,1)-V(0,0) ) );
      SLocal(0,1)= invvol * (  ( V(1,0)-V(1,1) )* V(1,1)          - ( V(0,1)-V(0,0) )*V(0,1) );
      SLocal(0,2)= invvol * ( -( V(1,0)-V(1,1) )* V(1,0)          + ( V(0,1)-V(0,0) )*V(0,0) );

      SLocal(1,1)=  invvol * (  V(1,1)*V(1,1) + V(0,1)*V(0,1) );
      SLocal(1,2)=  invvol * ( -V(1,1)*V(1,0) - V(0,1)*V(0,0) );

      SLocal(2,2)=  invvol * ( V(1,0)*V(1,0)+ V(0,0)*V(0,0) );

      SLocal(1,0)=SLocal(0,1);
      SLocal(2,0)=SLocal(0,2);
      SLocal(2,1)=SLocal(1,2);
      
      // Quadrature for heat transfer coefficient
      double KLocal=(kappa(cells(icell,0))+kappa(cells(icell,1))+kappa(cells(icell,2)))/3.0;

      // Assemble into global stiffness matrix
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
          SGlobal(cells(icell,i),cells(icell,j))+=KLocal*SLocal(i,j);
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

        SGlobal(i0,i0)+=h*alpha(ireg)/3.0;
        SGlobal(i0,i1)+=h*alpha(ireg)/6.0;

        SGlobal(i1,i0)+=h*alpha(ireg)/6.0;
        SGlobal(i1,i1)+=h*alpha(ireg)/3.0;
      }
        
    }
    // Return global stiffness matrix
    return pSGlobal;
  }

  // Assemble zero right hand side for mixed Dirichlet/Homogeneus Neumann problem
  std::shared_ptr<numcxx::TArray1<double>>
  assemble_general_heat_rhs(
    numcxx::SimpleGrid &grid,  // Discretization grid
    numcxx::DArray1& f,        // heat source (per node)
    numcxx::DArray1& g,        // boundary ambient temperature
    numcxx::DArray1& alpha     // boundary heat transfer coefficient (large value marks Dirichlet)
    )
  {
    auto ndim=grid.spacedim();
    int npoints=grid.npoints();
    int ncells=grid.ncells();
    int nbfaces=grid.nbfaces();
    auto cells=grid.get_cells();   // Local-global dof map
    auto points=grid.get_points();  
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();


    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim, ndim);
    auto &V=*pV;

    auto pRhs=numcxx::DArray1::create(npoints);
    auto &Rhs=*pRhs;

    Rhs=0.0;
    for (int icell=0; icell<ncells; icell++)
    {
      
      // Fill matrix V
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);

      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);


      // Compute determinant
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double fac=0.5*det/3.0;
      int i;
      i=cells(icell,0);
      Rhs(i)+= fac*f(i);

      i=cells(icell,1);
      Rhs(i)+= fac*f(i);

      i=cells(icell,2);
      Rhs(i)+= fac*f(i);

    }
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int ireg=bfaceregions(ibface);
      if (alpha(ireg)>=DirichletPenalty)
      {
        int i;

        i=bfaces(ibface,0);
        Rhs(i)+=DirichletPenalty*g(ireg);

        i=bfaces(ibface,1);
        Rhs(i)+=DirichletPenalty*g(ireg);
      }
      else
      {
        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        // first oder quadrature 
        int i;
        i=bfaces(ibface,0);
        Rhs(i)+=0.5*h*g(ireg);

        i=bfaces(ibface,1);
        Rhs(i)+=0.5*h*g(ireg);
      }
    }
    return pRhs;
  }
  
/////////////////////////////////////////////////////////
  void assemble_transient_heat_matrix_and_rhs(
    numcxx::SimpleGrid &grid,// Discretization grid
    numcxx::DSparseMatrix& S,
    numcxx::DArray1& Rhs,
    numcxx::DArray1& OldSol,
    numcxx::DArray1& f,    // heat source (per node)
    numcxx::DArray1& g,    // boundary ambient temperature
    numcxx::DArray1& kappa, // heat conduction coefficient (per node)
    numcxx::DArray1& alpha, // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    double tau, // time step
    double theta //  choice of method
    )
  {
    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    int npoints=grid.npoints();
    int ncells=grid.ncells();
    double tauinv=1.0/tau;
    
    // Local stiffness matrix
    auto pSLocal=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;
    
    // Local mass matrix
    auto pMLocal0=numcxx::DMatrix::create({{2,1,1},{1,2,1},{1,1,2}}); // from Stroud quadrature...
//    auto pMLocal0=numcxx::DMatrix::create({{4,0,0},{0,4,0},{0,0,4}}); // mass lumping
    auto &MLocal0=*pMLocal0;
    MLocal0*=1.0/12.0;
    
    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim, ndim);
    auto &V=*pV;
    
    // Loop over all elements (cells) of the triangulation
    Rhs=0.0;
    S(0,0)=0;
    S.clear();
    for (int icell=0; icell<ncells; icell++)
    {
      // Fill matrix V
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);
      
      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);
      

      // Compute determinant
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double vol = 0.5*det;
      double invvol=1.0/vol;
      
      // Compute entries of local stiffness matrix
      SLocal(0,0)= invvol * (  ( V(1,0)-V(1,1) )*( V(1,0)-V(1,1) )+( V(0,1)-V(0,0) )*( V(0,1)-V(0,0) ) );
      SLocal(0,1)= invvol * (  ( V(1,0)-V(1,1) )* V(1,1)          - ( V(0,1)-V(0,0) )*V(0,1) );
      SLocal(0,2)= invvol * ( -( V(1,0)-V(1,1) )* V(1,0)          + ( V(0,1)-V(0,0) )*V(0,0) );
      
      SLocal(1,1)=  invvol * (  V(1,1)*V(1,1) + V(0,1)*V(0,1) );
      SLocal(1,2)=  invvol * ( -V(1,1)*V(1,0) - V(0,1)*V(0,0) );
      
      SLocal(2,2)=  invvol * ( V(1,0)*V(1,0)+ V(0,0)*V(0,0) );
      
      SLocal(1,0)=SLocal(0,1);
      SLocal(2,0)=SLocal(0,2);
      SLocal(2,1)=SLocal(1,2);
      
      // Quadrature for heat transfer coefficient
      double KLocal=(kappa(cells(icell,0))+kappa(cells(icell,1))+kappa(cells(icell,2)))/3.0;
      
      
      double fac=vol/3.0;
      int i;
      for (int i=0;i<=ndim;i++)
      {
        Rhs(cells(icell,i))+= fac*f(cells(icell,i));
        for (int j=0;j<=ndim;j++)
        {
          S(cells(icell,i),cells(icell,j))+=theta*KLocal*SLocal(i,j)+ tauinv*vol*MLocal0(i,j);
          Rhs(cells(icell,j))+=((1.0-theta)*KLocal*SLocal(i,j)+ tauinv*vol*MLocal0(i,j))* OldSol(cells(icell,j));
        }
      }    
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
        S(i,i)+=alpha(ireg);
        Rhs(i)+=DirichletPenalty*g(ireg);
        
        i=bfaces(ibface,1);
        S(i,i)+=alpha(ireg);
        Rhs(i)+=DirichletPenalty*g(ireg);
        
      }
      else if (alpha(ireg)>0.0)
      {

        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        // first oder quadrature, needs to be replaced by second order
        int i0=bfaces(ibface,0);
        int i1=bfaces(ibface,1);
        
        S(i0,i0)+=theta*h*alpha(ireg)/3.0;
        S(i0,i1)+=theta*h*alpha(ireg)/6.0;
        
        S(i1,i0)+=theta*h*alpha(ireg)/6.0;
        S(i1,i1)+=theta*h*alpha(ireg)/3.0;
        
        Rhs(i0)+=OldSol(i0)*(1.0-theta)*h*alpha(ireg)/3.0;
        Rhs(i0)+=OldSol(i1)*(1.0-theta)*h*alpha(ireg)/6.0;
        Rhs(i1)+=OldSol(i0)*(1.0-theta)*h*alpha(ireg)/6.0;
        Rhs(i1)+=OldSol(i1)*(1.0-theta)*h*alpha(ireg)/3.0;
        
        Rhs(i0)+=0.5*h*g(ireg);
        Rhs(i0)+=0.5*h*g(ireg);
      }
    }
    S.flush();
  }

////////////////////////////////////////////////////////

  
  // Interfacing with python requires to work with smart pointers... 
    inline 
      std::shared_ptr<numcxx::DSparseMatrix> 
      assemble_heat_matrix(
        std::shared_ptr<numcxx::SimpleGrid> pg, 
        std::shared_ptr<numcxx::IArray1> diri
        )
    {return assemble_heat_matrix(*pg,*diri);}
    
    inline
  std::shared_ptr<numcxx::TArray1<double>>
  assemble_heat_rhs(std::shared_ptr<numcxx::SimpleGrid> pg, 
                    std::shared_ptr<numcxx::IArray1> diri,  
                    std::shared_ptr<numcxx::DArray1> bcval
    )
  {
    return assemble_heat_rhs(*pg,*diri,*bcval);
  }


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
  assemble_general_heat_rhs(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> f,        // heat source (per node)
    std::shared_ptr<numcxx::DArray1> g,        // boundary ambient temperature
    std::shared_ptr<numcxx::DArray1> alpha     // boundary heat transfer coefficient (large value marks Dirichlet)
    )
  {
    return assemble_general_heat_rhs(*grid,*f,*g,*alpha);
  }



  void assemble_transient_heat_matrix_and_rhs(
    std::shared_ptr<numcxx::SimpleGrid> grid,// Discretization grid
    std::shared_ptr<numcxx::DSparseMatrix> S,
    std::shared_ptr<numcxx::DArray1> Rhs,
    std::shared_ptr<numcxx::DArray1> OldSol,
    std::shared_ptr<numcxx::DArray1> f,    // heat source (per node)
    std::shared_ptr<numcxx::DArray1> g,    // boundary ambient temperature
    std::shared_ptr<numcxx::DArray1> kappa, // heat conduction coefficient (per node)
    std::shared_ptr<numcxx::DArray1> alpha, // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    double tau,
    double theta
    )
  {
    assemble_transient_heat_matrix_and_rhs(*grid,*S,*Rhs,*OldSol,*f,*g,*kappa,*alpha,tau,theta);
  }


}
#endif
