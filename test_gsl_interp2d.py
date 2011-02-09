from wrap_gsl_interp2d import Interp2DPeriodic

import numpy as np

from nose.tools import set_trace, ok_, eq_

def test_periodic():
    n, m = 10, 7
    nn = 512
    max_iv = float(nn)
    Y = np.linspace(0, max_iv, nn, endpoint=False)
    # print Y
    X = Y[:,np.newaxis]
    arr = np.sin(2. * np.pi * n * X / max_iv) + np.sin(2. * np.pi * m * Y / max_iv)

    # ensure the array edges are equal, for gsl periodic cspline purposes
    # ok_(np.allclose(arr[:,0], arr[:,-1]))
    # ok_(np.allclose(arr[0,:], arr[-1,:]))

    interp = Interp2DPeriodic(max_iv, max_iv, arr)

    ok_(np.allclose(arr[0,0], interp.eval(0.0, 0.0)))
    np.allclose(arr[0,0], interp.eval(110.0, 0.0))

    # tot_err = 0.0
    # for yval in range(len(Y)):
        # interps = [interp.eval(xx, Y[yval]) for xx in Y]
        # ok_(np.allclose(arr[:,yval], interps))
        # tot_err += np.sum(np.abs(arr[:,yval]-interps))
    # print tot_err

    XX1 = np.linspace(0, max_iv, nn, endpoint=False)
    interps1 = [interp.eval(xx, Y[3]) for xx in XX1]
    # XX2 = np.linspace(0, max_iv, 4*nn, endpoint=False)
    # interps2 = [interp.eval(xx, Y[3]) for xx in XX2]
    import pylab as pl
    pl.ion()
    # pl.imshow(arr)
    pl.plot(XX1, interps1, 'go')
    # pl.plot(XX2, interps2, 'rD')
    pl.plot(Y, arr[:, 3], 'bx-')
    set_trace()
