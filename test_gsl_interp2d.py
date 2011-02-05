from wrap_gsl_interp2d import Interp2DPeriodic

import numpy as np

from nose.tools import set_trace, ok_, eq_

def test_periodic():
    n, m = 10, 7
    nn = 100
    Y = np.linspace(0, nn+1, nn+2)
    print Y
    X = Y[:,np.newaxis]
    arr = np.sin(2. * np.pi * n * X / (nn+1)) + np.sin(2. * np.pi * m * Y / (nn+1))

    # ensure the array edges are equal, for gsl periodic cspline purposes
    ok_(np.allclose(arr[:,0], arr[:,-1]))
    ok_(np.allclose(arr[0,:], arr[-1,:]))

    interp = Interp2DPeriodic(Y, Y, arr)

    ok_(np.allclose(arr[0,0], interp.eval(0.0, 0.0)))
    np.allclose(arr[0,0], interp.eval(110.0, 0.0))

    tot_err = 0.0
    for yval in range(nn+1):
        interps = [interp.eval(xx, Y[yval]) for xx in Y]
        ok_(np.allclose(arr[:,yval], interps))
        tot_err += np.sum(np.abs(arr[:,yval]-interps))
    print tot_err
