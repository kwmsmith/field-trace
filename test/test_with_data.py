import tables
from sys import getsizeof

from wrap_gsl_interp2d import Interp2DPeriodic

import field_trace
import contour_tree as ct

import numpy as np
from itertools import izip

from test_critical_point_network import visualize

import pylab as pl
pl.ion()

def nth_timeslice(h5file, gpname, n):
    if isinstance(h5file, tables.file.File):
        dta = h5file
    elif isinstance(h5file, basestring):
        dta = tables.openFile(h5file, mode='r')
    gp = dta.getNode('/%s' % gpname)
    children = gp._v_children.keys()
    nth_ch = sorted(children)[n]
    return gp._v_children[nth_ch]

def h5_gen(h5file, gpname):
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

def test_region_contains():
    for psi_arr in h5_gen('data.h5', 'psi'):
        psi_arr = psi_arr.read()

def num_nonredundant_nulls(nulls, shape):
    from scipy.ndimage import label
    arr = np.zeros(shape, dtype=np.int32)
    for null in nulls:
        arr[null.loc] = 1
    larr, nlabels = label(arr)
    return nlabels

def save_figs():
    from upsample import upsample
    for bx, by, psi_arr in izip(h5_gen('data.h5', 'bx'),
                                h5_gen('data.h5', 'by'),
                                h5_gen('data.h5', 'psi')):
        bx = upsample(bx.read(), factor=1)
        by = upsample(by.read(), factor=1)
        psi_arr = upsample(psi_arr.read(), factor=1)
        # min_regions = field_trace.detect_min_regions(psi_arr, min_size=20)
        # mask = field_trace.regions_to_mask(psi_arr.shape, min_regions)
        # all_nulls, all_regions = field_trace.nulls_and_regions(psi_arr, chatty=True)
        all_nulls = field_trace.get_nulls(psi_arr)
        mins = [null for null in all_nulls if null.is_minimum()]
        maxs = [null for null in all_nulls if null.is_maximum()]
        saddles = [null for null in all_nulls if null.is_saddle()]
        print "-"*80
        print "num mins: %d unique mins: %d" % (len(mins), num_nonredundant_nulls(mins, psi_arr.shape))
        print "num maxs: %d unique maxs: %d" % (len(maxs), num_nonredundant_nulls(maxs, psi_arr.shape))
        print "num saddles: %d unique saddles: %d" % (len(saddles), num_nonredundant_nulls(saddles, psi_arr.shape))
        print "num nulls: %d unique nulls: %d" % (len(all_nulls), num_nonredundant_nulls(all_nulls, psi_arr.shape))
        # print "peaks - saddles: %d" % (len(maxs)+len(mins) - len(saddles))
        # print "num maxs: %d" % (len(maxs))
        # print "num mins: %d" % (len(mins))
        # print "num saddles: %d" % (len(saddles))
        # field_trace.find_region_ncontained(psi_arr.shape, all_regions)
        # n2regions = field_trace.regions_by_n_contained(all_regions)
        # mask = field_trace.regions_to_mask(psi_arr.shape, n2regions[0])
        # mask += field_trace.regions_to_mask(psi_arr.shape, n2regions[1])
        # mask += field_trace.regions_to_mask(psi_arr.shape, n2regions[2])
        # modb = np.sqrt(bx**2 + by**2)
        # print "saving to file"
        # field_trace.save_fig(modb, 'bmag_%03d' % ctr)
        # field_trace.save_fig(mask, 'mask_%03d' % ctr)
        # overlay = modb
        # overlay[mask] = overlay.max()
        # peak_x = [peak.x0 for peak in peaks]
        # peak_y = [peak.y0 for peak in peaks]
        # field_trace.save_fig_with_scatter(overlay, (peak_x, peak_y), 'bmagmaskpeaks_%03d' % ctr)
        # field_trace.save_fig(overlay, 'bmagmask_%03d' % ctr)
        # ctr += 1
        break

# save_figs()

def total_graph_memory(gr):
    tot_mem = 0
    try:
        tot_mem += getsizeof(gr.adj)
        tot_mem += getsizeof(gr.node)
        tot_mem += getsizeof(gr.edge)
    except AttributeError:
        pass
    return tot_mem

def test_contour_tree():
    from time import ctime
    for n in (0, -1):
        psi_arr = nth_timeslice('data.h5', 'psi', n=n)
        # for bx, by, psi_arr in izip(h5_gen('data.h5', 'bx'),
                                    # h5_gen('data.h5', 'by'),
                                    # h5_gen('data.h5', 'psi')):
        # bx = bx.read()
        # by = by.read()
        arr = psi_arr.read()
        print "array memory size: %d" % arr.nbytes
        def height_func(n):
            return (arr[n], n)
        print ctime(), "meshing array...",
        mesh = ct.make_mesh(arr)
        print ctime(), "done"
        print "mesh memory size: %d" % total_graph_memory(mesh)
        # num_eql = 0
        # for node in mesh:
            # for nbr in mesh.neighbors(node):
                # if arr[node] == arr[nbr]:
                    # num_eql += 1
        # print "number of equal height nodes: %d" % num_eql

        print ctime(), "ct.sparse_contour_tree()...",
        c_tree_sparse, regions_sparse = ct.sparse_contour_tree(mesh, height_func)
        print ctime(), "done"
        cpts_sparse = ct.critical_points(c_tree_sparse)
        peaks_sparse = cpts_sparse['peaks']
        passes_sparse = cpts_sparse['passes']
        pits_sparse = cpts_sparse['pits']
        print ctime(), "computing contour tree...",
        c_tree = ct.contour_tree(mesh, height_func)
        print ctime(), "done"
        print ctime(), "computing regions...",
        regions = ct.get_regions_full(c_tree)
        print ctime(), "done"
        print ctime(), "computing critical points...",
        cpts = ct.critical_points(c_tree)
        print ctime(), "done"
        peaks = cpts['peaks']
        passes = cpts['passes']
        pits = cpts['pits']
        print "sparse: peaks + pits - passes = %d" % (len(peaks_sparse) + len(pits_sparse) - len(passes_sparse))
        print "len(crit_pts_sparse) = %d" % (len(peaks_sparse) + len(pits_sparse) + len(passes_sparse))
        print "peaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        print "len(crit_pts) = %d" % (len(peaks) + len(pits) + len(passes))
        print "c_tree memory size: %d" % total_graph_memory(c_tree)
        print "getsizeof(regions): %d" % getsizeof(regions)
        regions2area = [len(r) for r in regions.values()]
        if 1:
            import pylab as pl
            pl.figure()
            pl.hist(regions2area, bins=pl.sqrt(len(regions2area)))
            raw_input("enter to continue")
            pl.close('all')
            filled_arr = arr.copy()
            def hf((a,b)): return height_func(a), height_func(b)
            ctr = 0
            for region in sorted(regions, key=hf, reverse=True):
                X = [_[0] for _ in regions[region]]
                Y = [_[1] for _ in regions[region]]
                filled_arr[X,Y] = 2*arr.max()
                if not ctr:
                    visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='gray', new_fig=False)
                    raw_input("enter to continue")
                ctr += 1
                ctr %= len(regions) / 20
            visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='gray', new_fig=False)
            raw_input("enter to continue")
            pl.close('all')

        del c_tree
        del regions
        del cpts
        del c_tree_sparse
        del regions_sparse

test_contour_tree()

"""
n=-1:
    (array([858, 192,  95,  59,  50,  28,  22,  12,  16,   9,   4,   3,   4,
             1,   8,   3,   2,   4,   4,   0,   0,   0,   2,   0,   1,   2,
             0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   2]),
     array([  2.00000000e+00,   1.86756757e+02,   3.71513514e+02,
             5.56270270e+02,   7.41027027e+02,   9.25783784e+02,
             1.11054054e+03,   1.29529730e+03,   1.48005405e+03,
             1.66481081e+03,   1.84956757e+03,   2.03432432e+03,
             2.21908108e+03,   2.40383784e+03,   2.58859459e+03,
             2.77335135e+03,   2.95810811e+03,   3.14286486e+03,
             3.32762162e+03,   3.51237838e+03,   3.69713514e+03,
             3.88189189e+03,   4.06664865e+03,   4.25140541e+03,
             4.43616216e+03,   4.62091892e+03,   4.80567568e+03,
             4.99043243e+03,   5.17518919e+03,   5.35994595e+03,
             5.54470270e+03,   5.72945946e+03,   5.91421622e+03,
             6.09897297e+03,   6.28372973e+03,   6.46848649e+03,
             6.65324324e+03,   6.83800000e+03]))
"""
