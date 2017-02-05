#!/usr/bin/env python
title="Simple test for P1 finite elements with heat source"

import numcxx
from numcxx import numcxxplot
import numpy, math
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
        [0.7,0.7],
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

bcfac=numcxx.asdarray([fem2d.Dirichlet,fem2d.Dirichlet,fem2d.Dirichlet,fem2d.Dirichlet,fem2d.Dirichlet])

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


do_plot=False
norms=[]
def solve(h):
    vol=h*h*math.sqrt(3.0)/4.0
    geom.set_regionvolumes([vol])

    grid=numcxx.SimpleGrid.create(geom,"zpaAqDQ")

    nnodes=grid.npoints()
    points=grid.get_points()

    source=numcxx.DArray1.create(nnodes)
    exact=numcxx.DArray1.create(nnodes)
    kappa=numcxx.DArray1.create(nnodes)
    solution=numcxx.DArray1.create(nnodes)
    error=numcxx.DArray1.create(nnodes)

    for i in range(nnodes):
        kappa[i]=1.0
        x=points[i][0]
        y=points[i][1]
        exact[i]=math.sin(math.pi*x)*math.sin(math.pi*y)
        source[i]=2.0*math.pi*math.pi*exact[i]

    S=numcxx.DSparseMatrix.create(nnodes,nnodes);
    Solver=numcxx.DSolverUMFPACK.create(S)
    Rhs=numcxx.DArray1.create(nnodes);
    fem2d.assemble_heat_problem_with_source(grid,bcfac,bcval,source,kappa,S, Rhs);

    

    Solver.update()
    Solver.solve(solution,Rhs)

    for i in range(nnodes):
        error[i]=exact[i]-solution[i]
    l2norm=fem2d.l2norm(grid,error)
    h1norm=fem2d.h1norm(grid,error)
    thisnorms=[nnodes,h,l2norm,h1norm,numcxx.norm2(error)/math.sqrt(nnodes), numcxx.normi(error)]
    print(thisnorms)
    norms.append(thisnorms)

    if do_plot:
        toplot=error
        triang=numcxxplot.triangulation(grid)
        plt.tricontourf(triang, numcxx.asnumpy(toplot),20,cmap='gnuplot')
        plt.title(title)
        plt.colorbar()
        plt.tricontour(triang, numcxx.asnumpy(toplot),20,colors="black",linestyles="solid",colorbar='none')
        plt.show()




for h in [ 0.1*2**(-i) for i in range(6)]:
     solve(h)


plt.figure(1)
plt.clf()
plt.grid()
norms=numpy.array(norms)
oh=norms[:,1]
ohh=oh**2

plt.loglog(norms[:,1],norms[:,2],'r-o',label='$L^2$')
plt.loglog(norms[:,1],norms[:,3],'g-o',label='$H^1$ (seminorm)')
plt.loglog(norms[:,1],norms[:,4],'b-o',label='discr. $L^2$')
plt.loglog(norms[:,1],norms[:,5],'m-o',label='$L^\infty$')
plt.loglog(norms[:,1],0.25*oh,'g:',linewidth=3,label='$0.25h$')
plt.loglog(norms[:,1],0.25*ohh,'m:',linewidth=3,label='$0.25h^2$')

plt.ylabel('error norm')
plt.xlabel('h')

plt.legend(loc='upper left')
plt.draw()
plt.savefig("test-fem2d-convergence-h.pdf")


plt.clf()
plt.grid()
norms=numpy.array(norms)
oh=norms[:,1]
ohh=oh**2

plt.loglog(norms[:,0],norms[:,2],'r-o',label='$L^2$')
plt.loglog(norms[:,0],norms[:,3],'g-o',label='$H^1$ (seminorm)')
plt.loglog(norms[:,0],norms[:,4],'b-o',label='discr. $L^2$')
plt.loglog(norms[:,0],norms[:,5],'m-o',label='$L^\infty$')
plt.loglog(norms[:,0],0.25*oh,'g:',linewidth=3,label='$0.25h$')
plt.loglog(norms[:,0],0.25*ohh,'m:',linewidth=3,label='$0.25h^2$')

plt.ylabel('Error norm')
plt.xlabel('Number of unknowns')


plt.legend(loc='lower left')
plt.draw()
plt.savefig("test-fem2d-convergence-N.pdf")

