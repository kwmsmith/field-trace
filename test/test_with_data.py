from wrap_gsl_interp2d import Interp2DPeriodic

import field_trace

import numpy as np


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

psi_arrs = [arr.read() for arr in h5_gen('data.h5', 'psi')]
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

    nulls = field_trace.find_and_cull_cells(dd)
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
    # N = 64
    # n, m = 3, 5
    # test_data = tv.sin_cos_arr(N, n, m)
    test_data = psi_arrs[-1].astype(np.double)
    N = test_data.shape[0]

    dd = field_trace.Derivator(test_data, N, N)
    nulls = field_trace.find_and_cull_cells(dd)
    saddles = [null for null in nulls if null.is_saddle()]
    peaks = [null for null in nulls if not null.is_saddle()]
    saddle0s = [(s.x0, s.y0) for s in saddles]
    peak0s = [(p.x0, p.y0) for p in peaks]

    x0, y0 = saddle0s[0]

    psi_interp = Interp2DPeriodic(N, N, test_data)
    level_val = psi_interp.eval(x0, y0)

    print "computing level sets"
    null2level = field_trace.level_sets(test_data, psi_interp, saddles[:50])

    if 1:
        import pylab as pl
        pl.ion()
        all_masks = field_trace.marked_to_mask(test_data.shape, null2level.values())
        masked_data = test_data.copy()
        masked_data[all_masks] = test_data.max()
        pl.imshow(masked_data, cmap='hot', interpolation='nearest')
        X, Y = zip(*saddle0s)
        pl.scatter(Y, X, c='k')
        X, Y = zip(*peak0s)
        pl.scatter(Y, X, c='b')
        raw_input("enter to continue")
        # for level in null2level.values():
            # masked_data = test_data.copy()
            # mask = field_trace.marked_to_mask(test_data.shape, [level])
            # masked_data[mask] = test_data.max()
            # pl.imshow(masked_data, cmap='hot', interpolation='nearest')
            # X, Y = zip(*saddle0s)
            # pl.scatter(Y, X, c='k')
            # X, Y = zip(*peak0s)
            # pl.scatter(Y, X, c='b')
            # raw_input("enter to continue")

test_level_set()
