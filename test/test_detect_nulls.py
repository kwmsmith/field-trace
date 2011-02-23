import numpy as np
from kaw_analysis import test_vcalc as tv
from kaw_analysis import vcalc

from wrap_gsl_interp2d import Interp2DPeriodic

from nose.tools import ok_, eq_, set_trace

import field_trace

from scipy import ndimage

global_N = 64

def _test_find_null_2x2():

    class _dummy(object): pass
    dd = _dummy()

    oeigsystem = field_trace.eigsystem
    try:
        field_trace.eigsystem = lambda x,y,z: (None, None)

        dd.perp_deriv1 = np.array([[-2.,-2.],[1.,1.]])
        dd.perp_deriv2 = np.array([[-1., 2.],[-1., 2]])
        null = field_trace.find_null_cells(dd)[0]
        ok_(np.allclose((null.x0, null.y0), (2./3, 1./3)))

        dd.perp_deriv1 = np.array([[-1.,-1.],[1.,1.]])
        dd.perp_deriv2 = np.array([[-1., 1.],[-1., 1.]])
        null = field_trace.find_null_cells(dd)[0]
        ok_(np.allclose((null.x0, null.y0), (1./2, 1./2)))
    finally:
        field_trace.eigsystem = oeigsystem

def test_find_nulls():
    def tester(n, m):
        N = global_N
        test_data = tv.sin_cos_arr(N, n, m)
        dd = field_trace.Derivator(test_data, N, N)
        nulls = field_trace.find_null_cells(dd)
        eq_(4*n*m, len(nulls))
        cell_locs = [n.loc for n in nulls]
        null_locs = [(n.x0, n.y0) for n in nulls]

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
        for m in range(3, 5):
            yield tester, n, m

def test_derivator():
    N = global_N
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

def test_level_set():
    N = 64
    n, m = 1, 4
    test_data = tv.sin_cos_arr(N, n, m)

    dd = field_trace.Derivator(test_data, N, N)
    nulls = field_trace.find_null_cells(dd)
    saddles = [null for null in nulls if null.is_saddle()]
    peaks = [null for null in nulls if not null.is_saddle()]
    saddle0s = [(s.x0, s.y0) for s in saddles]
    peak0s = [(p.x0, p.y0) for p in peaks]

    x0, y0 = saddle0s[0]

    psi_interp = Interp2DPeriodic(N, N, test_data)
    level_val = psi_interp.eval(x0, y0)

    level_sets = [null.levelset for null in nulls]

    mask = field_trace.marked_to_mask(test_data.shape, level_sets)

    masked_data = test_data.copy()
    masked_data[mask] = test_data.max()

    if 0:
        import pylab as pl
        pl.ion()
        pl.imshow(masked_data, cmap='hot', interpolation='nearest')
        # # Plot the grid points
        # if 0:
            # X = np.linspace(0, dta.shape[0]-1, dta.shape[0])
            # for i in range(dta.shape[0]):
                # Y = np.zeros(dta.shape[1])
                # Y.fill(i)
                # pl.scatter(Y, X, c='m')
        X, Y = zip(*saddle0s)
        pl.scatter(Y, X, c='k')
        X, Y = zip(*peak0s)
        pl.scatter(Y, X, c='b')
        raw_input("enter to continue")

def test_partition_regions():
    boundary = np.array([(i, 3) for i in range(10)] + [(i, 7) for i in range(10)])
    boundary = field_trace.Region(boundary[:,0], boundary[:,1])
    nx, ny = 10, 10
    regions = field_trace.partition_regions(nx, ny, boundary)
    eq_(len(regions), 2)

    arr = np.zeros((nx, ny), dtype=np.int_)
    field_trace.set_region(arr, boundary, 0)
    field_trace.set_region(arr, regions[0], 1)
    field_trace.set_region(arr, regions[1], 2)

    for row in arr:
        ok_(np.allclose(row, [1, 1, 1, 0, 2, 2, 2, 0, 1, 1]), `row`)

    if 0:
        import pylab as pl
        pl.ion()
        pl.imshow(arr, cmap='hot', interpolation='nearest')
        raw_input("enter to continue")


def test_neighbors4():
    eq_(field_trace.neighbors4(0, 0, 10, 10), ((1, 0), (0, 1), (9, 0), (0, 9)))

class test_cc_label(object):

    def setup(self):
        self.PLOT = False
        self.background = 0
        self.empty = 1
        self.square_hole = np.zeros((10, 10), dtype=np.int32)
        self.square_hole.fill(self.empty)
        self.square_hole[3:8,3:8] = self.background
        self.square_hole[4:7,4:7] = self.empty

    def test_cc_label(self):
        larr = field_trace.connected_component_label(self.square_hole)

        ok_(np.all(larr == \
                 np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 0, 0, 0, 0, 0, 1, 1],
                           [1, 1, 1, 0, 2, 2, 2, 0, 1, 1],
                           [1, 1, 1, 0, 2, 2, 2, 0, 1, 1],
                           [1, 1, 1, 0, 2, 2, 2, 0, 1, 1],
                           [1, 1, 1, 0, 0, 0, 0, 0, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=np.int32)))

        if 0:
            import pylab as pl
            pl.ion()
            pl.imshow(larr, cmap='hot', interpolation='nearest')
            raw_input("enter to continue")

    def test_square_hole_shifted(self):
        self.square_hole = np.roll(self.square_hole, shift=4, axis=1)
        larr = field_trace.connected_component_label(self.square_hole)

        ok_(np.all(larr == \
                 np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                           [2, 0, 1, 1, 1, 1, 1, 0, 2, 2],
                           [2, 0, 1, 1, 1, 1, 1, 0, 2, 2],
                           [2, 0, 1, 1, 1, 1, 1, 0, 2, 2],
                           [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=np.int32)))

        if 0:
            import pylab as pl
            pl.ion()
            pl.close('all')
            pl.imshow(self.square_hole, interpolation='nearest')
            pl.figure()
            pl.imshow(larr, interpolation='nearest')
            raw_input("enter to continue")

    def test_square_hole_shifted_up(self):
        self.square_hole = np.roll(self.square_hole, shift=4, axis=0)
        outarr = self.square_hole.copy()
        larr = field_trace.connected_component_label(self.square_hole, output=outarr)

        eq_(id(outarr), id(larr))

        ok_(np.all(larr == \
                 np.array([[1, 1, 1, 0, 2, 2, 2, 0, 1, 1],
                           [1, 1, 1, 0, 0, 0, 0, 0, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 0, 0, 0, 0, 0, 1, 1],
                           [1, 1, 1, 0, 2, 2, 2, 0, 1, 1],
                           [1, 1, 1, 0, 2, 2, 2, 0, 1, 1]], dtype=np.int32)))

        if 0:
            import pylab as pl
            pl.ion()
            pl.close('all')
            pl.imshow(self.square_hole, interpolation='nearest')
            pl.figure()
            pl.imshow(larr, interpolation='nearest')
            raw_input("enter to continue")
