#!/usr/bin/env python
title="Demonstrate transient heat source"
import numcxx
from numcxx import numcxxplot
import numpy,math
import fem2d
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

bcfac=numcxx.asdarray([0,fem2d.DirichletPenalty,fem2d.DirichletPenalty,fem2d.DirichletPenalty,fem2d.DirichletPenalty])

bcval=numcxx.asdarray([0,0,0,0,0])

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
triang=numcxxplot.triangulation(grid)


T=1.0
N=200
tau=1.0/N
theta=0.5


nnodes=grid.npoints()
points=grid.get_points()

Matrix=numcxx.DSparseMatrix.create(nnodes,nnodes);
Solver=numcxx.DSolverUMFPACK.create(Matrix)

source=numcxx.DArray1.create(nnodes)
kappa=numcxx.DArray1.create(nnodes)
Rhs=numcxx.DArray1.create(nnodes)
Sol=numcxx.DArray1.create(nnodes)
OldSol=numcxx.DArray1.create(nnodes)

for i in range(nnodes):
    source[i]=0.0
    kappa[i]=1.0
    Sol[i]=0.0
    
sigma=1.0e-2
t=0.0
for n in range(10*N):
    t=t+tau;
    xpos=0.5+0.4*numpy.cos(20*t)
    ypos=0.5+0.4*numpy.sin(20*t)
    print("time step %d"%n)
    for i in range(nnodes):
        dx=points[i][0]-xpos
        dy=points[i][1]-ypos
        source[i]=0.01*math.exp((dx*dx+dy*dy)/sigma*sigma)
        OldSol[i]=Sol[i]

    fem2d.assemble_transient_heat_matrix_and_rhs(
        grid,
        Matrix,
        Rhs,
        OldSol,
        source,
        bcval,
        kappa,
        bcfac,
        tau,
        theta)

    Solver.update()
    Solver.solve(Sol,Rhs)
    npsol=numcxx.asnumpy(Sol)
    plt.clf()
    plt.tricontourf(triang, npsol,20,cmap='gnuplot')
    plt.title("%s, t=%f\n"%(title,t))
    plt.colorbar()
    plt.tricontour(triang, npsol,20,colors="black",linestyles="solid",colorbar='none')
    plt.pause(1.0e-10)

plt.show()
