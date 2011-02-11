from wrap_gsl_interp2d import Interp2DPeriodic

import numpy as np

from nose.tools import set_trace, ok_, eq_

global_NN = 128

def make_arr(NN, modes_x, modes_y, max_xy):
    Y = np.linspace(0, max_xy, NN, endpoint=False)
    X = Y[:,np.newaxis]
    arr = (np.sin(2. * np.pi * modes_x * X / max_xy) +
            np.sin(2. * np.pi * modes_y * Y / max_xy))
    return (arr, X, Y)


def test_periodic():

    # wavenumbers in X & Y dimensions, resp.
    n, m = 10, 3 
    max_iv = 3.1415926

    arr, X, Y = make_arr(global_NN, n, m, max_iv)

    interp = Interp2DPeriodic(max_iv, max_iv, arr)

    ok_(np.allclose(arr[0,0], interp.eval(0.0, 0.0)),
            "%f != %f" % (arr[0,0], interp.eval(0.0, 0.0)))

    tot_err = 0.0
    for yval in range(len(Y)):
        interps = [interp.eval(xx, Y[yval]) for xx in Y]
        ok_(np.allclose(arr[:,yval], interps))
        tot_err += np.sum(np.abs(arr[:,yval]-interps))
    ok_(np.allclose(0.0, tot_err))

    # test refinement...
    for yval in range(len(Y)):
        # sample at twice the array resolution.
        XX_h = np.linspace(0, max_iv, 2*global_NN, endpoint=False)
        interp_h = [interp.eval(xx, Y[yval]) for xx in XX_h]
        # sample at 4X the array res.
        XX_h2 = np.linspace(0, max_iv, 4*global_NN, endpoint=False)
        interp_h2 = [interp.eval(xx, Y[yval]) for xx in XX_h2]
        # compare the values where they (should) overlap
        ok_(np.allclose(interp_h, interp_h2[::2]))

    if 0:
        PLOT_Y = 3
        XX1 = np.linspace(0, max_iv, global_NN, endpoint=False)
        interps1 = [interp.eval(xx, Y[PLOT_Y]) for xx in XX1]
        # interps_y = [interp.eval(15.0, yy) for yy in Y]
        XX2 = np.linspace(0, max_iv, 4*global_NN, endpoint=False)
        interps2 = [interp.eval(xx, Y[PLOT_Y]) for xx in XX2]
        import pylab as pl
        pl.ion()
        # pl.imshow(arr)
        # pl.plot(Y, interps_y, 'go')
        pl.plot(XX2, interps2, 'rD')
        pl.plot(XX1, interps1, 'ko')
        pl.plot(Y, arr[:, PLOT_Y], 'bx-')
        raw_input("--> enter to quit")

def test_many_interps():
    max_xy = 10.0
    arr, X, Y = make_arr(global_NN, modes_x=10, modes_y=20, max_xy=max_xy)

    ninterps = 3
    interpolators = [Interp2DPeriodic(max_xy, max_xy, arr) for _ in range(ninterps)]

    tot_err = 0.0
    for yval in range(len(Y)):
        vecs = [np.array([interp.eval(xx, Y[yval]) for xx in Y]) for interp in interpolators]
        for vec in vecs:
            ok_(np.allclose(vecs[0], vec))
            tot_err += np.sum(np.abs(vecs[0]-vec))
    ok_(np.allclose(0.0, tot_err))
