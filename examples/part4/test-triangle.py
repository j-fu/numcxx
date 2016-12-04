#!/usr/bin/env python
import numcxx
import numpy
import matplotlib
from matplotlib import pyplot as plt



geom=numcxx.Geometry.create()
geom.set_points([[0,0.1],[1,0],[1.1,1]])
geom.set_bfaces([[0,1],[1,2],[2,0]])
geom.set_bfaceregions([1,2,3])
geom.set_regionpoints([[0.5,0.25]])
geom.set_regionvolumes([0.001])
geom.set_regionnumbers([1])

print(numcxx.asnumpy(geom.bfaceregions))



numcxx.plot.plotGeometry(plt,geom)
plt.show()
grid=numcxx.SimpleGrid.create(geom,"zpaAq")










numcxx.plot.plotGrid(plt,grid)
plt.show()
print(matplotlib.get_backend())
