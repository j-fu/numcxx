#ifndef FEM2D_H
#define FEM2D_H

#include <functional>

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

  void  assemble_transient_heat_problem(
    const numcxx::SimpleGrid &grid,
    const numcxx::DArray1& bcfac,
    const numcxx::DArray1& bcval,
    const numcxx::DArray1& source,
    const numcxx::DArray1& kappa,
    double tau, // time step size
    numcxx::DArray1 &OldSol,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Rhs);

  
  void  assemble_and_apply_nonlinear_heat(
    const numcxx::SimpleGrid &grid,
    const numcxx::DArray1& bcfac,
    const numcxx::DArray1& bcval,
    const numcxx::DArray1& source,
    std::function <void(const double, double&, double&)> fkappa,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Sol,
    numcxx::DArray1 &Rhs);



  void initialize_bc(
    numcxx::SimpleGrid &grid,// Discretization grid
    numcxx::DArray1& g,
    numcxx::DArray1& Sol
    );

  double l2norm(const numcxx::SimpleGrid &grid, 
                const numcxx::DArray1 &u);
  
  double h1norm(const numcxx::SimpleGrid &grid, 
                const numcxx::DArray1 &u);
}
#endif


