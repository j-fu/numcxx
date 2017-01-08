#ifndef FEM2D_H
#define FEM2D_H

#include <numcxx/simplegrid.hxx>
#include <numcxx/tsparsematrix.hxx>
#include <cmath>
#include <iostream>
namespace fem2d
{
  inline std::shared_ptr<numcxx::TSparseMatrix<double>> 
  assemble_heat_matrix(
    numcxx::SimpleGrid &g,
    numcxx::IArray1 &bcdirichlet
    )
  {
    auto ndim=g.spacedim();
    auto points=g.get_points();
    auto cells=g.get_cells();
    int npoints=g.npoints();
    int ncells=g.ncells();

    
    auto pSLocal=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto pSGlobal=numcxx::DSparseMatrix::create(npoints, npoints);
    auto pV=numcxx::DMatrix::create(ndim, ndim);

    auto &SLocal=*pSLocal;
    auto &V=*pV;
    auto &SGlobal=*pSGlobal;
//    std::cout << points << std::endl;
//    std::cout << cells << std::endl;

    for (int icell=0; icell<ncells; icell++)
    {
      
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);

      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);

//      std::cout << V << std::endl;
      
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double invdet = 1.0/det;

      
      SLocal(0,0)= invdet * (  ( V(1,0)-V(1,1) )*( V(1,0)-V(1,1) )+( V(0,1)-V(0,0) )*( V(0,1)-V(0,0) ) );
      SLocal(0,1)= invdet * (  ( V(1,0)-V(1,1) )* V(1,1)          - ( V(0,1)-V(0,0) )*V(0,1) );
      SLocal(0,2)= invdet * ( -( V(1,0)-V(1,1) )* V(1,0)          + ( V(0,1)-V(0,0) )*V(0,0) );

      SLocal(1,1)=  invdet * (  V(1,1)*V(1,1) + V(0,1)*V(0,1) );
      SLocal(1,2)=  invdet * ( -V(1,1)*V(1,0) - V(0,1)*V(0,0) );

      SLocal(2,2)=  invdet * ( V(1,0)*V(1,0)+ V(0,0)*V(0,0) );

      SLocal(1,0)=SLocal(0,1);
      SLocal(2,0)=SLocal(0,2);
      SLocal(2,1)=SLocal(1,2);

//      std::cout << SLocal << std::endl;
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
          SGlobal(cells(icell,i),cells(icell,j))+=SLocal(i,j);
    }    


    int nbfaces=g.nbfaces();
    auto bfaces=g.get_bfaces();
    auto bfaceregions=g.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int ireg=bfaceregions(ibface);
      if (bcdirichlet(ireg))
      {
        int i;
        i=bfaces(ibface,0);
        SGlobal(i,i)+=1.0e30;
        i=bfaces(ibface,1);
        SGlobal(i,i)+=1.0e30;
      }
    }
    return pSGlobal;
  }

  std::shared_ptr<numcxx::TArray1<double>>
  assemble_heat_rhs(
    numcxx::SimpleGrid &g,
    numcxx::IArray1 &bcdirichlet,
    numcxx::DArray1 &bcval
    )
  {
    int npoints=g.npoints();

    auto pRhs=numcxx::DArray1::create(npoints);
    auto &Rhs=*pRhs;
    Rhs=0.0;
    int nbfaces=g.nbfaces();
    auto bfaces=g.get_bfaces();
    auto bfaceregions=g.get_bfaceregions();

    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int ireg=bfaceregions(ibface);
      if (bcdirichlet(ireg))
      {
        int i;
        i=bfaces(ibface,0);
        Rhs(i)+=1.0e30*bcval(ireg);
        i=bfaces(ibface,1);
        Rhs(i)+=1.0e30*bcval(ireg);
      }
    }
    return pRhs;
  }
  
  
  // Interfacing with python requires to work with smart pointers... 
  inline std::shared_ptr<numcxx::DSparseMatrix> 
  assemble_heat_matrix(std::shared_ptr<numcxx::SimpleGrid> pg, std::shared_ptr<numcxx::IArray1> diri)
  {return assemble_heat_matrix(*pg,*diri);}
  
  inline
  std::shared_ptr<numcxx::TArray1<double>>
  assemble_heat_rhs(std::shared_ptr<numcxx::SimpleGrid> pg, 
                    std::shared_ptr<numcxx::IArray1> diri,  
                    std::shared_ptr<numcxx::DArray1> bcval)
  {
    return assemble_heat_rhs(*pg,*diri,*bcval);
  }


}
#endif
