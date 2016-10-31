#!/usr/bin/env python

import numcxx
import numpy



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

xa=numcxx.asnumcxx(x)

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

