#!/usr/bin/env python
#  Python script to test triangle mesh generator via numcxx
import numcxx
from numcxx import numcxxplot
import numpy
import mymodule


# Create a geometry object
geom=numcxx.Geometry.create()

# Give initial list of points
geom.set_points(
    [
        [0,0.1],
        [1,0],
        [1.1,1],
    ])

# Give initial list of boundary faces.
# Indices refer to position in point list
geom.set_bfaces(
    [
        [0,1],
        [1,2],
        [2,0]
    ])

# Give region numbers for boundary faces
# These should be larger than 0
geom.set_bfaceregions([1,2,3])

# Give interior region points to mark
# region
geom.set_regionpoints([[0.5,0.25]])

# Give maximal area of triangles
# in region corresponding to region volumes
geom.set_regionvolumes([0.1])

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




f=mymodule.myfunction(grid)
nf=numcxx.asnumpy(f)
print(nf)

