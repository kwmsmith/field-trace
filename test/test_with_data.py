from kaw_analysis import test_vcalc as tv
from wrap_gsl_interp2d import Interp2DPeriodic

import field_trace

import numpy as np
from itertools import izip

import pylab as pl
pl.ion()


def h5_gen(h5file, gpname):
    import sys
    import tables
    if isinstance(h5file, tables.file.File):
        dta = h5file
    elif isinstance(h5file, basestring):
        dta = tables.openFile(h5file, mode='r')
    gp = dta.getNode('/%s' % gpname)
    for arr in gp:
        yield arr
    dta.close()

# psi_arrs = [arr.read() for arr in h5_gen('data.h5', 'psi')]
# cur_arrs = [arr.read() for arr in h5_gen('data.h5', 'cur')]

def upsample(arr, fac):
    new_a = np.empty((arr.shape[0] * fac, arr.shape[1] * fac), dtype=np.double)

    Xs = np.linspace(0.0, arr.shape[0], new_a.shape[0], endpoint=False)
    Ys = np.linspace(0.0, arr.shape[1], new_a.shape[1], endpoint=False)

    interp = Interp2DPeriodic(arr.shape[0], arr.shape[1], arr.astype(np.double))
    
    for idx, X in enumerate(Xs):
        new_a[idx, :] = [interp.eval(X, Y) for Y in Ys]

    return new_a

def test_tracer():
    from kaw_analysis import test_vcalc as tv
    N = 64
    dta = tv.sin_cos_prod(N, 5, 5)
    # dta = psi_arrs[-1]

    dd = field_trace.Derivator(dta, dta.shape[0], dta.shape[1])

    nulls = field_trace.find_null_cells(dd)
    saddles = [null for null in nulls if null.is_saddle()]
    peaks = [null for null in nulls if not null.is_saddle()]

    error = 0.0
    for null in nulls:
        error += abs(dd.perp_deriv1_interp.eval(null.x0, null.y0))
        error += abs(dd.perp_deriv2_interp.eval(null.x0, null.y0))
    error /= len(nulls)

    saddle0s = [(s.x0, s.y0) for s in saddles]

    peak0s = [(p.x0, p.y0) for p in peaks]

    # finished, traces = field_trace.make_skeleton(
            # arr=dta, deriv=dd, saddles=saddles, ncomb=1, start_offset=1.e-3, hit_radius=1.e-3)

    traces = []
    for saddle in saddles:
        for outgoing in (True, False):
            for sgn in (1,-1):
                trace = field_trace.trace_from_saddle(
                        arr=dta, deriv=dd, start_saddle=saddle,
                        start_offset=sgn*1.e-2, outgoing=outgoing, timeout=0.001)
                traces.append(trace)

    if 1:
        import pylab as pl
        pl.ion()
        pl.imshow(dta, cmap='hot')
        for tr in traces:
            tr = np.array(tr)
            for i in range(0, len(tr[0]), 2):
                X,Y = tr[:,i], tr[:,i+1]
                pl.scatter(Y, X, c='r')

        # Plot the grid points
        if 0:
            X = np.linspace(0, dta.shape[0]-1, dta.shape[0])
            for i in range(dta.shape[0]):
                Y = np.zeros(dta.shape[1])
                Y.fill(i)
                pl.scatter(Y, X, c='m')
            
        X, Y = zip(*saddle0s)
        pl.scatter(Y, X, c='k')
        X, Y = zip(*peak0s)
        pl.scatter(Y, X, c='b')
        import pdb; pdb.set_trace()

# test_tracer()

def test_level_set():
    # N = 128
    # n, m = 10, 15
    # test_data = tv.sin_cos_arr(N, n, m)
    arr_idx = -2
    test_data = psi_arrs[arr_idx].astype(np.double)
    N = test_data.shape[0]

    print "locating nulls"
    dd = field_trace.Derivator(test_data, N, N)
    nulls = field_trace.find_null_cells(dd)
    null0s = [(n.x0, n.y0) for n in nulls]

    psi_interp = Interp2DPeriodic(N, N, test_data)

    print "computing level sets"
    levels = [null.levelset for null in nulls]

    print "classifying nulls"

    peaks = []
    saddles = []
    for null in nulls:
        if null.is_saddle():
            saddles.append(null)
        else:
            peaks.append(null)

    peak0s = [(p.x0, p.y0) for p in peaks]
    saddle0s = [(s.x0, s.y0) for s in saddles]

    print "getting min regions"
    regions = []
    for null in nulls:
        regions.extend(null.regions)

    print "number of regions: %d" % len(regions)

    min_regions = field_trace.filter_min_regions(regions)

    print "number of min regions: %d" % len(min_regions)

    print "plotting"

    if 1:
        import pylab as pl
        pl.ion()
        pl.figure()
        dta = np.zeros(test_data.shape, dtype=np.int32)
        for min_region in min_regions:
            dta[min_region.xs, min_region.ys] = 1
        pl.imshow(dta, interpolation='nearest')
        pl.figure()
        pl.imshow(cur_arrs[arr_idx])
        pl.figure()
        pl.imshow(psi_arrs[arr_idx])
        raw_input("enter to continue")

        
    if 0:
        import pylab as pl
        pl.ion()
        all_masks = field_trace.marked_to_mask(test_data.shape, levels)
        masked_data = test_data.copy()
        masked_data[all_masks] = test_data.max()
        pl.imshow(masked_data, cmap='hot', interpolation='nearest')
        X, Y = zip(*peak0s)
        pl.scatter(Y, X, c='k')
        X, Y = zip(*saddle0s)
        pl.scatter(Y, X, c='b')

        if 0:
            for level in levels:
                masked_data = test_data.copy()
                mask = field_trace.marked_to_mask(test_data.shape, [level])
                masked_data[mask] = test_data.max()
                pl.imshow(masked_data, cmap='hot', interpolation='nearest')
                X, Y = zip(*null0s)
                pl.scatter(Y, X, c='k')
                X, Y = zip(*null0s)
                pl.scatter(Y, X, c='b')

        raw_input("enter to continue")

def test_detect_min_regions():
    from itertools import izip

    for bx, by, psi_arr in izip(h5_gen('data.h5', 'bx'),
                                h5_gen('data.h5', 'by'),
                                h5_gen('data.h5', 'psi')):
        bx = bx.read()
        by = by.read()
        psi_arr = psi_arr.read()
        min_regions = field_trace.detect_min_regions(psi_arr, min_size=20)
        mask = field_trace.regions_to_mask(psi_arr.shape, min_regions)
        pl.figure()
        pl.imshow(bx**2 + by**2, cmap='hot')
        pl.figure()
        pl.imshow(mask, cmap='hot', interpolation='nearest')
        raw_input("enter to continue")
        pl.close('all')

def save_figs():
    ctr = 0
    for bx, by, psi_arr in izip(h5_gen('data.h5', 'bx'),
                                h5_gen('data.h5', 'by'),
                                h5_gen('data.h5', 'psi')):
        bx = bx.read()
        by = by.read()
        psi_arr = psi_arr.read()
        min_regions = field_trace.detect_min_regions(psi_arr, min_size=20)
        mask = field_trace.regions_to_mask(psi_arr.shape, min_regions)
        print "saving to file"
        field_trace.save_fig(bx**2 + by**2, 'bmag_%03d' % ctr)
        field_trace.save_fig(mask, 'mask_%03d' % ctr)
        ctr += 1

save_figs()
