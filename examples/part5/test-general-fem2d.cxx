#include <numcxx/numcxx.hxx>
#include "fem2d.hxx"


int main(void)
{
    auto geom=numcxx::Geometry::create();
    geom->set_points(
        {
            {0,0},
            {1,0},
            {1,1},
            {0,1.5},
        });
    
    geom->set_bfaces(
        {
            {0,1},
            {1,2},
            {2,3},
            {3,0}
        });
    geom->set_bfaceregions({1,2,3,4});

    auto bcfac=numcxx::DArray1::create({0,fem2d::Dirichlet,0,fem2d::Dirichlet,0});
    auto bcval=numcxx::DArray1::create({0,10,0,0,0});

    geom->set_regionpoints({{0.5,0.55}});
    geom->set_regionvolumes({0.1});
    geom->set_regionnumbers({1});

    auto grid=numcxx::SimpleGrid::create(*geom,"zpaAqD");

    auto nnodes=grid->npoints();

    auto source=numcxx::DArray1::create(nnodes);
    auto kappa=numcxx::DArray1::create(nnodes);

    *source=0.0;
    *kappa=1.0;

    auto S=numcxx::DSparseMatrix::create(nnodes,nnodes);
    auto Rhs=numcxx::DArray1::create(nnodes);
    
    fem2d::assemble_heat_problem_with_source(grid,bcfac,bcval,source,kappa,S, Rhs);
    auto Solver=numcxx::DSolverUMFPACK::create(S);

    Solver->update();
    auto Sol=Rhs->copy();
    Solver->solve(Sol,Rhs);

    std::cout << *Sol << std::endl;
    
}

