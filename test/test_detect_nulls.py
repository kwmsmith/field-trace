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

def test_hessian_esys():
    def tester(n, m):
        N = global_N
        test_data = tv.sin_cos_prod(N, n, m)
        # from upsample import upsample
        # test_data = upsample(test_data, factor=4)
        dd = field_trace.Derivator(test_data)
        arr1 = np.zeros((N, N), dtype=np.double)
        arr2 = np.zeros((N, N), dtype=np.double)
        for i in range(N):
            for j in range(N):
                evals, evecs = field_trace.hessian_esys(dd, i, j)
                arr1[i,j] = min(evals)
                arr2[i,j] = max(evals)
        nulls = field_trace.find_null_cells_minimize(dd)
        null_locs = [(null.x0, null.y0) for null in nulls]
        saddles = [_ for _ in nulls if _.is_saddle()]
        maxs = [_ for _ in nulls if _.is_maximum()]
        mins = [_ for _ in nulls if _.is_minimum()]
        saddleXs, saddleYs = [_.x0 for _ in saddles], [_.y0 for _ in saddles]
        minXs, minYs = [_.x0 for _ in mins], [_.y0 for _ in mins]
        maxXs, maxYs = [_.x0 for _ in maxs], [_.y0 for _ in maxs]

        if 0:
            def ims(x):
                pl.imshow(x, interpolation='nearest', cmap='hot')
                pl.scatter(saddleYs, saddleXs, marker='x', c='k')
                pl.scatter(minYs, minXs, marker='o', c='b')
                pl.scatter(maxYs, maxXs, marker='d', c='r')
            import pylab as pl
            pl.ion()
            # ims(dd.deriv11)
            # pl.figure()
            # ims(dd.deriv22)
            # pl.figure()
            # ims(dd.deriv12)
            # pl.figure()
            ims(test_data)
            pl.figure()
            ims(arr1)
            pl.figure()
            ims(arr2)
            pl.figure()
            ims(arr1*arr2)
            raw_input("enter to continue")
            pl.close('all')

        # eq_(4*n*m, len(nulls))

    for n in (3, 7, 11):
        for m in (5, 9, 10):
            yield tester, n, m

def test_max_min_saddle():
    def tester(n, m):
        N = global_N
        test_data = tv.sin_cos_prod(N, n, m)
        from upsample import upsample
        test_data = upsample(test_data, factor=4)
        dd = field_trace.Derivator(test_data)
        nulls = field_trace.find_null_cells_minimize(dd)
        null_locs = [(null.x0, null.y0) for null in nulls]
        saddles = [_ for _ in nulls if _.is_saddle()]
        maxs = [_ for _ in nulls if _.is_maximum()]
        mins = [_ for _ in nulls if _.is_minimum()]
        eq_(len(saddles + maxs + mins), len(nulls))

        print ("4*n*m: %d num_saddles: %d num_max+num_min: %d" %
                (4*n*m, len(saddles), len(maxs+mins)))

        if 1:
            import pylab as pl
            pl.ion()
            pl.imshow(test_data)
            X = [_.x0 for _ in saddles]
            Y = [_.y0 for _ in saddles]
            pl.scatter(Y, X, marker='x', c='k')
            X = [_.x0 for _ in maxs]
            Y = [_.y0 for _ in maxs]
            pl.scatter(Y, X, c='r', marker='d')
            X = [_.x0 for _ in mins]
            Y = [_.y0 for _ in mins]
            pl.scatter(Y, X, c='b', marker='o')
            raw_input("enter to continue")
            pl.close('all')

    for n in (3, 7, 11):
        for m in (5, 9, 10):
            yield tester, n, m


def test_find_nulls():
    def tester(n, m):
        N = global_N
        test_data = tv.sin_cos_arr(N, n, m)
        from upsample import upsample
        test_data = upsample(test_data, factor=4)
        dd = field_trace.Derivator(test_data)
        nulls = field_trace.find_null_cells_minimize(dd)
        cell_locs = [null.loc for null in nulls]
        null_locs = [(null.x0, null.y0) for null in nulls]
        print 4*n*m, len(nulls)

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

        # eq_(4*n*m, len(nulls))

    for n in (3, 7, 11):
        for m in (5, 9, 10):
            yield tester, n, m

def test_derivator():
    N = global_N
    n, m = 10, 3
    test_data = tv.sin_cos_arr(N, n, m)
    dd = field_trace.Derivator(test_data)

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

    dd = field_trace.Derivator(test_data)
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
