import numpy as np
from kaw_analysis import test_vcalc as tv
from kaw_analysis import vcalc

from wrap_gsl_interp2d import Interp2DPeriodic

from nose.tools import ok_, eq_, set_trace

import field_trace

def test_find_null_2x2():
    class _dummy(object): pass
    dd = _dummy()
    dd.perp_deriv1 = np.array([[-2.,-2.],[1.,1.]])
    dd.perp_deriv2 = np.array([[-1., 2.],[-1., 2]])
    nulls = field_trace.find_null_cells(dd, step=1)
    ok_(np.allclose(nulls[0].find_null(), (2./3, 1./3)))

    dd.perp_deriv1 = np.array([[-1.,-1.],[1.,1.]])
    dd.perp_deriv2 = np.array([[-1., 1.],[-1., 1.]])
    nulls = field_trace.find_null_cells(dd, step=1)
    ok_(np.allclose(nulls[0].find_null(), (1./2, 1./2)))

def test_find_nulls():
    def tester(n, m):
        N = 256
        test_data = tv.sin_cos_arr(N, n, m)
        dd = field_trace.Derivator(test_data, N, N)
        nulls = field_trace.find_and_cull_cells(dd)
        eq_(4*n*m, len(nulls))
        cell_locs = [n.loc for n in nulls]
        null_locs = [n.find_null() for n in nulls]

        if 0:
            X, Y = zip(*null_locs)
            import pylab as pl
            pl.ion()
            pl.imshow(test_data)
            pl.scatter(Y, X)
            X, Y = zip(*cell_locs)
            pl.scatter(Y, X, c='r')
            raw_input("enter to continue")
            pl.close('all')

    for n in range(3, 5):
        for m in range(10, 12):
            yield tester, n, m

def _test_null_cell():
    nc = field_trace.NullCell((10, 11), 2, (512, 512))
    eq_(nc.covered(), set([(10,11),(11,11),(11,12),(10,12)]))
    nc = field_trace.NullCell((10, 11), 1, (512, 512))
    eq_(nc.covered(), set([(10,11)]))
    nc = field_trace.NullCell((10, 10), 2, (11,11))
    eq_(nc.covered(), set([(10,0),(0,10),(10,10),(0,0)]))

def test_eig_system():
    N = 256
    n, m = 3, 4
    test_data = tv.sin_cos_arr(N, n, m)

    dd = field_trace.Derivator(test_data, N, N)

    nulls = field_trace.find_and_cull_cells(dd)
    null_locs = [n.find_null() for n in nulls]

    for nl in null_locs:
        evals, evecs = field_trace.eigsystem(dd, nl[0], nl[1])

    for null in nulls:
        null.is_saddle()

    if 0:
        import pylab as pl
        pl.ion()
        pl.imshow(test_data)
        pl.figure()
        pl.imshow(psi_11)
        pl.title('XX')
        pl.figure()
        pl.imshow(psi_22)
        pl.title('YY')
        pl.figure()
        pl.imshow(psi_12)
        pl.title('XY')
        raw_input("enter to continue")

def test_derivator():
    N = 256
    n, m = 10, 3
    test_data = tv.sin_cos_arr(N, n, m)
    dd = field_trace.Derivator(test_data, N, N)

    psi_1 = vcalc.cderivative(test_data, 'X_DIR')
    psi_2 = vcalc.cderivative(test_data, 'Y_DIR')

    ok_(np.allclose(dd.deriv1, psi_1))
    ok_(np.allclose(dd.deriv2, psi_2))

    ok_(np.allclose(dd.perp_deriv1, psi_2))
    ok_(np.allclose(dd.perp_deriv2, -psi_1))

    psi_11 = vcalc.cderivative(test_data, 'X_DIR', order=2)
    psi_22 = vcalc.cderivative(test_data, 'Y_DIR', order=2)
    psi_12 = vcalc.cderivative(test_data, 'X_DIR')
    psi_12 = vcalc.cderivative(psi_12,    'Y_DIR')

    Xlen = float(N)

    psi_11_interp = Interp2DPeriodic(Xlen, Xlen, psi_11)
    psi_22_interp = Interp2DPeriodic(Xlen, Xlen, psi_22)
    psi_12_interp = Interp2DPeriodic(Xlen, Xlen, psi_12)

    num_pts = 100
    rand_Xs = np.random.uniform(0.0, float(N), num_pts)
    rand_Ys = np.random.uniform(0.0, float(N), num_pts)

    def compare_interps(dd_interp, interp, xys):
        dd_interps = [dd_interp.eval(x, y) for (x, y) in zip(*xys)]
        interps    =    [interp.eval(x, y) for (x, y) in zip(*xys)]
        ok_(np.allclose(dd_interps, interps))

    compare_interps(dd.deriv12_interp, psi_12_interp, [rand_Xs, rand_Ys])
    compare_interps(dd.deriv11_interp, psi_11_interp, [rand_Xs, rand_Ys])
    compare_interps(dd.deriv22_interp, psi_22_interp, [rand_Xs, rand_Ys])
