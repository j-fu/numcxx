#!/usr/bin/env python
import numcxx
import numpy

n=4

print("numcxx array to numpy:")
a=numcxx.DArray1.create(n)
for i in range(n):
    a[i]=i
npa=numcxx.asnumpy(a)
print(npa)


print("numpy array as numcxx darray:")
b=numpy.ndarray(n)
for i in range(n):
    b[i]=n-i

ncb=numcxx.asdarray(b)
for i in range(n):
    print(ncb[i])

print("numpy array as numcxx iarray:")
c=numpy.ndarray(n,dtype=numpy.intc)
for i in range(n):
    c[i]=n-i

ncc=numcxx.asiarray(c)
for i in range(n):
    print(ncc[i])


# Create a numcxx matrix from python list
M0=numcxx.asdmatrix([[2,3,5],
       [3,2,3],
       [9,5,7.0]]);

# Create right hand side from python list
f=numcxx.asdarray([1,2,3])
x=numcxx.DArray1.create(3)

# Create LU solver
LU=numcxx.DSolverLapackLU.create(M0)
# Solve
LU.solve(x,f)

# Print result
print(numcxx.asnumpy(x))


