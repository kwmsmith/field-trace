import numpy as np
from kaw_analysis import test_vcalc as tv
from kaw_analysis import vcalc

from nose.tools import ok_, eq_, set_trace

import field_trace

def test_find_null_2x2():
    ax = np.array([[-2.,-2.],[1.,1.]])
    ay = np.array([[-1., 2.],[-1., 2]])
    nulls = field_trace.find_null_cells(ax, ay, step=1)
    np.allclose(nulls[0].find_null(), (2./3, 1./3))

    ax = np.array([[-1.,-1.],[1.,1.]])
    ay = np.array([[-1., 1.],[-1., 1.]])
    nulls = field_trace.find_null_cells(ax, ay, step=1)
    np.allclose(nulls[0].find_null(), (1./2, 1./2))

def test_find_nulls():
    def tester(n, m):
        test_data = tv.sin_cos_arr(256, n, m)
        d_x = vcalc.cderivative(test_data, 'X_DIR')
        d_y = vcalc.cderivative(test_data, 'Y_DIR')
        nulls = field_trace.find_and_cull_cells(d_x, d_y)
        eq_(4*n*m, len(nulls))
        cell_locs = [n.loc for n in nulls]
        null_locs = [n.find_null() for n in nulls]
        # X, Y = zip(*null_locs)
        # import pylab as pl
        # pl.ion()
        # pl.imshow(test_data)
        # pl.scatter(Y, X)
        # X, Y = zip(*cell_locs)
        # pl.scatter(Y, X, c='r')
        # set_trace()

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
