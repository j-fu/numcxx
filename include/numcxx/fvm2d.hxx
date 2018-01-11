#ifndef FEM2D_H
#define FEM2D_H

#include <numcxx/simplegrid.hxx>
#include <numcxx/tsparsematrix.hxx>

namespace fvm2d
{
  
  /// BC value marking Dirichlet boundary condition
  const double Dirichlet=1.0e30;
  

  void  assemble_heat_problem(
    const numcxx::SimpleGrid &Grid, // Discretization grid containing triangulation
    const numcxx::DArray1& BCfac,   // Array of boudary factors  (per boundary region)
    const numcxx::DArray1& BCval,   // Array of boudary values  (per boundary region)
    const numcxx::DArray1& Source,  // Array of source values (per node)
    const numcxx::DArray1& Kappa,   // Array of heat coefficient values (per node)
    numcxx::DSparseMatrix &SGlobal, // Global stiffness matrix
    numcxx::DArray1 &Rhs);          // Right hand side
  
  // double l2norm(const numcxx::SimpleGrid &grid, 
  //               const numcxx::DArray1 &u);
  
  // double h1norm(const numcxx::SimpleGrid &grid, 
  //               const numcxx::DArray1 &u);
}
#endif


