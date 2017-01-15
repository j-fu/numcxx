#!/usr/bin/env python
#  Python script to test triangle mesh generator via numcxx
import numcxx
from numcxx import numcxxplot
import numpy
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
        [0,1.5],
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

diri=numcxx.asiarray([0,1,0,1,0])

dirival=numcxx.asdarray([0,10,0,0,0])

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



S=fem2d.assemble_heat_matrix(grid,diri)

Rhs=fem2d.assemble_heat_rhs(grid,diri,dirival)


Solver=numcxx.DSolverUMFPACK.create(S)
Solver.update()
Sol=Rhs.copy()
Solver.solve(Sol,Rhs)

npsol=numcxx.asnumpy(Sol)

# numcxxplot.plotGrid(plt,grid)
# print("Close window to finish!");
# plt.show()


triang=numcxxplot.triangulation(grid)
plt.tricontourf(triang, npsol,20,cmap='gnuplot')
plt.colorbar()
plt.tricontour(triang, npsol,20,colors="black",linestyles="solid",colorbar='none')


plt.show()