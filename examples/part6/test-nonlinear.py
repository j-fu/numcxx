#!/usr/bin/env python
title="Test transient  problem"
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

bcfac=numcxx.asdarray([0,fvm2d.DirichletPenalty,0,fvm2d.DirichletPenalty,0])

bcval=numcxx.asdarray([0,10,0,0,0])

# Give interior region points to mark
# region
geom.set_regionpoints([[0.5,0.55]])

# Give maximal area of triangles
# in region corresponding to region volumes
geom.set_regionvolumes([0.01])

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
numcxxplot.plotGrid(plt,grid)
print("Close window to finish!");
plt.show()
triang=numcxxplot.triangulation(grid)



nnodes=grid.npoints()


Matrix=numcxx.DSparseMatrix.create(nnodes,nnodes);
Solver=numcxx.DSolverUMFPACK.create(Matrix)
Res=numcxx.DArray1.create(nnodes)
Upd=numcxx.DArray1.create(nnodes)
Sol=numcxx.DArray1.create(nnodes)


for i in range(nnodes):
    Sol[i]=0.0
    Res[i]=0.0

fvm2d.initialize_bc(grid,Sol,bcval)



iter=0
norm=1
while iter<100 and norm >1.0e-13:
    fvm2d.assemble_and_apply_nonlinear_heat(grid,Matrix, Sol, Res,bcfac,bcval)
    Solver.update()
    
    Solver.solve(Upd,Res)
    oldnorm=norm
    norm=numcxx.norm2(Upd)
    print("norm=%8.4e contract=%8.5e"%(norm,norm/oldnorm))
    
    for i in range(nnodes):
        Sol[i]=Sol[i]-Upd[i]
    iter=iter+1

plt.clf()
npsol=numcxx.asnumpy(Sol)
plt.tricontourf(triang, npsol,20,cmap='gnuplot')
plt.colorbar()
plt.tricontour(triang, npsol,20,colors="black",linestyles="solid",colorbar='none')
plt.pause(1.0e-10)
plt.show()

    


