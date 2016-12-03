#!/usr/bin/env python

import numcxx
import numpy
import matplotlib
from matplotlib import pyplot as plt



a=numcxx.DArray1.create(4)

a.itemset(3,5.0)

print(a.item(3))
print(a[3])
print(a.shape(0))


npa=numcxx.DArray1AsNumpy(a)
print(npa[3])


x=numpy.array([3,4,5.0])
print(x.shape[0])



xa=numcxx.NumpyAsDArray1(x)
print(xa[1])

xa=numcxx.asnumcxx_darray(x)
na=numcxx.asnumpy(xa)

print(na[1])


a2=numcxx.DArray2.create(4,5)

for i in range(4):
    for j in range(5):
        a2[i][j]=10*i+j

print(a2.shape(0),a2.shape(1))


a2x=numcxx.asnumpy(a2)
for i in range(4):
    for j in range(5):
        print("%d %d %8.2g"%(i,j,a2x[i][j]))


M0=numpy.array(
        [   [2,3,5],
            [3,2,3],
            [9,5,7.0]]);
print(M0)

M1=numcxx.asnumcxx_matrix(M0)
MI=M1.calculate_inverse();
M2=numcxx.asnumpy(MI)
print(M2)

X=[1,2,3]
XX=numcxx.asnumcxx_iarray(X)

geom=numcxx.Geometry.create()
geom.set_points([[0,0.1],[1,0],[1.1,1]])
geom.set_bfaces([[0,1],[1,2],[2,0]])
geom.set_bfaceregions([1,2,3])
geom.set_regionpoints([[0.5,0.25]])
geom.set_regionvolumes([0.01])
geom.set_regionnumbers([1])

print(numcxx.asnumpy(geom.bfaceregions))


def plot_geom(plt,geom):
    points=geom.get_points()
    bfaces=geom.get_bfaces()
    bfaceregions=geom.get_bfaceregions()
    regionpoints=geom.get_regionpoints()
    regionnumbers=geom.get_regionnumbers()

    colortable=[ 
        (0,0,0),
        (1,0,0),
        (0,1,0),
        (0,0,1),
        (1,1,0),
        (0,1,1),
        (1,0,1)]

    for ipoint in range(regionpoints.shape[0]):
        plt.plot([regionpoints[ipoint][0]],[regionpoints[ipoint][1]],
                 'o-',color=colortable[regionnumbers[ipoint]%len(colortable)])

    for ibface in range(bfaces.shape[0]):
        p=[bfaces[ibface][0],bfaces[ibface][1]]
        plt.plot([points[p[0]][0],points[p[1]][0]],
                 [points[p[0]][1],points[p[1]][1]],
                 'o-',color=colortable[bfaceregions[ibface]%len(colortable)])

plot_geom(plt,geom)
plt.waitforbuttonpress()
grid=numcxx.SimpleGrid.create(geom,"zpaAq")


def plot_grid(plt,grid):
    points=grid.get_points()
    cells=grid.get_cells()
    bfaces=grid.get_bfaces()
    bfaceregions=grid.get_bfaceregions()

    colortable=[ 
        (0,0,0),
        (1,0,0),
        (0,1,0),
        (0,0,1),
        (1,1,0),
        (0,1,1),
        (1,0,1)]

    for icell in range(grid.ncells()):
        p=[cells[icell][0],cells[icell][1],cells[icell][2]]
        plt.plot([points[p[0]][0],points[p[1]][0],points[p[2]][0],points[p[0]][0]],
                 [points[p[0]][1],points[p[1]][1],points[p[2]][1],points[p[0]][1]],
                 'go-')

    for ibface in range(grid.nbfaces()):
        p=[bfaces[ibface][0],bfaces[ibface][1]]
        plt.plot([points[p[0]][0],points[p[1]][0]],
                 [points[p[0]][1],points[p[1]][1]],
                 'o-',color=colortable[bfaceregions[ibface]%len(colortable)])



def triangulation(grid):
    xy=grid.get_points()
    x=xy[:,[0]].squeeze()
    y=xy[:,[1]].squeeze()
    triangles=grid.get_cells()
    return matplotlib.tri.Triangulation(x,y,triangles)






#plt.triplot(triangulation(grid), 'bo-')

plt.clf()
plot_grid(plt,grid)
plt.waitforbuttonpress()
