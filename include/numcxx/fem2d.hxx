#ifndef FEM2D_H
#define FEM2D_H

#include <numcxx/simplegrid.hxx>
#include <numcxx/tsparsematrix.hxx>

namespace fem2d
{
  
  const double Dirichlet=1.0e30;
  
  void  assemble_heat_problem(
    const numcxx::SimpleGrid &Grid,
    const numcxx::DArray1& BCfac,
    const numcxx::DArray1& BCval,
    const numcxx::DArray1& Source,
    const numcxx::DArray1& Kappa,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Rhs);
  
  double l2norm(const numcxx::SimpleGrid &grid, 
                const numcxx::DArray1 &u);
  
  double h1norm(const numcxx::SimpleGrid &grid, 
                const numcxx::DArray1 &u);
}
#endif


