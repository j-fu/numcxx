import numpy
from matplotlib import pyplot
X=numpy.loadtxt("20-plot-X.dat");
U=numpy.loadtxt("20-plot-U.dat");


pyplot.figure(1)

pyplot.title("test plot")
pyplot.ylabel('U')
pyplot.xlabel('X')
pyplot.grid(True)
pyplot.plot(X,U,'green')
pyplot.savefig("20-numcxx-matplotlib.png")
pyplot.show()
pyplot.waitforbuttonpress()

pyplot.pause(1.0e-10)
