import matplotlib
from matplotlib import pyplot
from math import *
import numpy
import numcxx
import convdiff


for N in [20,40,80]:

    D=0.01
    V=1

    X=numpy.linspace(0,1,N)

    Uexp=numcxx.asnumpy(convdiff.solve_expfit(N,V,D))
    Ucnt=numcxx.asnumpy(convdiff.solve_central(N,V,D))
    Uupw=numcxx.asnumpy(convdiff.solve_upwind(N,V,D))
    

    pyplot.figure(1)
    pyplot.clf()
    pyplot.grid()
    pyplot.ylim(-0.5,1)
    pyplot.plot(X,Uexp,label='expfit')
    pyplot.plot(X,Ucnt,label='central')
    pyplot.plot(X,Uupw,label='upwind')
    pyplot.legend(loc='upper left')
    pyplot.draw()
    pyplot.savefig("convdiff-%d.pdf"%N)







