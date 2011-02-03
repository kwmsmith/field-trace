from bilinear import bilinear
import numpy as np

from nose.tools import ok_, eq_, set_trace

def test_arr():
    I = np.linspace(0, 10, 200, endpoint=False)
    Y = np.sin(I)
    X = Y[:,np.newaxis]

    arr = X+Y

    for i in range(arr.shape[0]-1):
        for j in range(arr.shape[1]-1):
            f00 = arr[i,j]
            f10 = arr[i+1,j]
            f01 = arr[i,j+1]
            f11 = arr[i+1,j+1]
            eq_(arr[i,j], bilinear(f00, f10, f01, f11, 0.0, 0.0))
            eq_(arr[i+1,j], bilinear(f00, f10, f01, f11, 1.0, 0.0))
            eq_(arr[i,j+1], bilinear(f00, f10, f01, f11, 0.0, 1.0))
            eq_(arr[i+1,j+1], bilinear(f00, f10, f01, f11, 1.0, 1.0))

# import pylab as pl
# pl.imshow(arr)
# pl.show()
