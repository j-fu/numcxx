#!/usr/bin/env python
title="Simple test for 2D FVM implementation"
import numcxx
from numcxx import numcxxplot
import numpy
import fvm2d
import matplotlib
from matplotlib import pyplot as plt


# Create a geometry object
geom=numcxx.Geometry.create()

# Give initial list of points
geom.set_points(
    [
        [0,0],
        [1,0],
        [1,1],
        [0,1],
    ])

# Give initial list of boundary faces.
# Indices refer to position in point list
geom.set_bfaces(
    [
        [0,1],
        [1,2],
        [2,3],
        [3,0]
    ])



# Give region numbers for boundary faces
# These should be larger than 0
geom.set_bfaceregions([1,2,3,4])

bcfac=numcxx.asdarray([0,fvm2d.Dirichlet,0,fvm2d.Dirichlet,0])

bcval=numcxx.asdarray([0,10,0,0,0])

# Give interior region points to mark
# region
geom.set_regionpoints([[0.5,0.55]])

# Give maximal area of triangles
# in region corresponding to region volumes
geom.set_regionvolumes([0.001])

# Give region numbers for regions corresponding
# to region point
geom.set_regionnumbers([1])


# Create triangulation
# Flags:
# -z  Numbers all items starting from zero (rather than one).
# -p  Triangulates a Planar Straight Line Graph
# -a  Applies a maximum triangle area constraint.
# -A  Applies attributes to identify triangles in certain regions.
# -q  Quality mesh generation.  A minimum angle may be specified.
# -D  Conforming Delaunay:  all triangles are truly Delaunay.
grid=numcxx.SimpleGrid.create(geom,"zpaAqD")

nnodes=grid.npoints()
source=numcxx.DArray1.create(nnodes)
kappa=numcxx.DArray1.create(nnodes)
for i in range(nnodes):
    source[i]=0.0
    kappa[i]=1.0

nnodes=grid.npoints()


Rhs=numcxx.DArray1.create(nnodes)
Sol=numcxx.DArray1.create(nnodes)
Matrix=numcxx.DSparseMatrix.create(nnodes,nnodes);

fvm2d.assemble_heat_problem_with_source(grid,bcfac,bcval,source,kappa,Matrix, Rhs);


Solver=numcxx.DSolverUMFPACK.create(Matrix)
Solver.update()
Solver.solve(Sol,Rhs)



npsol=numcxx.asnumpy(Sol)



triang=numcxxplot.triangulation(grid)
plt.tricontourf(triang, npsol,20,cmap='plasma')
plt.title(title)
plt.colorbar()
plt.tricontour(triang, npsol,20,colors="black",linestyles="solid",colorbar='none')


plt.show()
