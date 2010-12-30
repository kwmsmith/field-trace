from wrap_gsl_interp import Interp, Interp2D

import numpy as np
# x[i] = i + 0.5 * sin (i);
# y[i] = i + cos (i * i);

I = np.arange(10)
X = I + 0.5 * np.sin(I)
Y = I + np.cos(I*I)

ii = Interp(X, Y)

XIs = np.linspace(0.0, 10.0, num=10.0/0.01)

YIs = [ii.eval(xi) for xi in XIs]

# for (x,y) in zip(XIs, YIs):
    # print "%g %g" % (x,y)

N = 400
arr = np.random.rand(N,N)

i2d = Interp2D(arr)

for i in range(len(arr)):
    for j in range(len(arr[0])):
        assert np.allclose(i2d.eval(i, j), arr[i,j]), `(i,j,i2d.eval(i,j),arr[i,j])`
