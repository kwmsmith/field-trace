from __future__ import print_function
import tables
from sys import getsizeof
import pickle

from wrap_gsl_interp2d import Interp2DPeriodic

import field_trace
import contour_tree as ct

import numpy as np
from itertools import izip

from test_critical_point_network import visualize

import pylab as pl

import _critical_points as _cp

from copy import copy

def nth_timeslice(h5file, gpname, n):
    if isinstance(h5file, tables.file.File):
        dta = h5file
    elif isinstance(h5file, basestring):
        dta = tables.openFile(h5file, mode='r')
    gp = dta.getNode('/%s' % gpname)
    children = gp._v_children.keys()
    nth_ch = sorted(children)[n]
    arr = gp._v_children[nth_ch].read()
    dta.close()
    return arr

def h5_gen(h5file, gpnames):
    if isinstance(gpnames, basestring):
        gpnames = (gpnames,)
    if isinstance(h5file, tables.file.File):
        dta = h5file
    elif isinstance(h5file, basestring):
        dta = tables.openFile(h5file, mode='r')
    gps = [dta.getNode('/%s' % gpname) for gpname in gpnames]
    for arrs in izip(*gps):
        aas = [arr.read() for arr in arrs]
        yield aas
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
    from kaw_analysis import test_vcalc as tv#{{{
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
        import pdb; pdb.set_trace()#}}}

# test_tracer()

def test_level_set():
    # N = 128#{{{
    # n, m = 10, 15
    # test_data = tv.sin_cos_arr(N, n, m)
    arr_idx = -2
    test_data = psi_arrs[arr_idx].astype(np.double)
    N = test_data.shape[0]

    print("locating nulls")
    dd = field_trace.Derivator(test_data, N, N)
    nulls = field_trace.find_null_cells(dd)
    null0s = [(n.x0, n.y0) for n in nulls]

    psi_interp = Interp2DPeriodic(N, N, test_data)

    print("computing level sets")
    levels = [null.levelset for null in nulls]

    print("classifying nulls")

    peaks = []
    saddles = []
    for null in nulls:
        if null.is_saddle():
            saddles.append(null)
        else:
            peaks.append(null)

    peak0s = [(p.x0, p.y0) for p in peaks]
    saddle0s = [(s.x0, s.y0) for s in saddles]

    print( "getting min regions")
    regions = []
    for null in nulls:
        regions.extend(null.regions)

    print( "number of regions: %d" % len(regions))

    min_regions = field_trace.filter_min_regions(regions)

    print( "number of min regions: %d" % len(min_regions))

    print( "plotting")

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

        raw_input("enter to continue")#}}}

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
    from upsample import upsample#{{{
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
        print( "-"*80)
        print( "num mins: %d unique mins: %d" % (len(mins), num_nonredundant_nulls(mins, psi_arr.shape)))
        print( "num maxs: %d unique maxs: %d" % (len(maxs), num_nonredundant_nulls(maxs, psi_arr.shape)))
        print( "num saddles: %d unique saddles: %d" % (len(saddles), num_nonredundant_nulls(saddles, psi_arr.shape)))
        print( "num nulls: %d unique nulls: %d" % (len(all_nulls), num_nonredundant_nulls(all_nulls, psi_arr.shape)))
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
        break#}}}

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

def logger(ss, newline=True):
    from time import ctime
    from sys import stderr
    stderr.write("%s: " % ctime())
    stderr.write(ss)
    if newline:
        stderr.write('\n')

def compare_cp_and_contour_tree():
    logger("loading arrays...")#{{{
    ctr = 0
    ppp = []
    for (psi_arr,) in h5_gen('data.h5', ('psi',)):
        # psi_arr = nth_timeslice('data.h5', 'psi', -1)
        # bx_arr =  nth_timeslice('data.h5', 'bx',  -1)
        # by_arr =  nth_timeslice('data.h5', 'by',  -1)
        logger("getting topological surface... %d" % ctr)
        surf = _cp.TopoSurface(psi_arr)
        scpts = surf.crit_pts
        logger("done")
        ppp.append((ctr, scpts.peaks, scpts.pits, scpts.passes))
        logger("%d + %d - %d = %d" % (len(scpts.peaks), len(scpts.pits), len(scpts.passes),
                                      len(scpts.peaks) + len(scpts.pits) - len(scpts.passes)))
        pickle.dump(ppp, open('topo_surface_cpts.dat', 'wb'))
        ctr += 1
        # def height_func(n):
            # return (psi_arr[n], n)
        # logger("meshing psi_arr...")
        # mesh = ct.make_mesh(psi_arr)
        # logger("getting contour tree...")
        # c_tree = ct.contour_tree(mesh, height_func)
        # logger("getting contour tree critical points...")
        # cpts = ct.critical_points(c_tree)
        # logger("done...")
        # import pdb; pdb.set_trace()#}}}

def test_contour_tree():
    ctr = 0#{{{
    flux_tube_areas_record = []
    for bx_arr, by_arr, psi_arr in h5_gen('data.h5', ('bx', 'by', 'psi')):
        b_mag = np.sqrt(bx_arr**2 + by_arr**2)
        logger("array memory size: %d" % psi_arr.nbytes)
        def height_func(n):
            return (psi_arr[n], n)
        logger("meshing array...")
        mesh = ct.make_mesh(psi_arr)
        logger("done")
        logger("mesh memory size: %d" % total_graph_memory(mesh))

        logger("computing contour tree...")
        c_tree = ct.contour_tree(mesh, height_func)
        logger("done")

        def region_func(r):
            return len(r)

        logger("pruning small regions...")
        ct.prune_regions(
                c_tree,
                region_func=region_func,
                threshold=3,
                )
        logger("done")

        logger("computing critical points...")
        cpts = ct.critical_points(c_tree)
        logger("done")

        logger("condensing small regions...")
        ct.prune_regions(
                c_tree,
                region_func=region_func,
                threshold=psi_arr.size/200,
                )
        logger("done")

        logger("computing regions...")
        regions = ct.get_regions(c_tree)
        logger("done")

        flux_tubes = ct.get_flux_tubes(c_tree)
        flux_tube_areas = np.array([len(ft) for (e, ft) in flux_tubes], dtype=np.int64)
        flux_tube_areas_record.append(flux_tube_areas)

        pickle.dump(flux_tube_areas_record, open('flux_tube_areas.dat', 'wb'))

        flux_tube_mask = np.zeros(psi_arr.shape, dtype=np.bool_)
        for (e, ft) in flux_tubes:
            for p in ft:
                flux_tube_mask[p] = True

        peaks = cpts.peaks
        passes = cpts.passes
        pits = cpts.pits
        all_cpts = peaks.union(passes).union(pits)

        logger("peaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes)))
        logger("len(crit_pts) = %d" % (len(peaks) + len(pits) + len(passes)))

        cpt_grads = [b_mag[cpt] / b_mag.min() for cpt in all_cpts]
        logger("b_mag.mean()=%f" % b_mag.mean())
        logger("b_mag.std()=%f" % b_mag.std())
        logger("b_mag.max()=%f" % b_mag.max())
        logger("b_mag.min()=%f" % b_mag.min())
        if 1:
            pl.figure()
            pl.hist(cpt_grads, bins=pl.sqrt(len(cpt_grads)))
            pl.title("cpoint gradient values")
            pl.figure()
            pl.hist(flux_tube_areas, bins=pl.sqrt(len(flux_tube_areas)))
            pl.title("flux tube areas")
            pl.figure()
            filled_arr = psi_arr.copy()
            filled_arr[flux_tube_mask] = 2*psi_arr.max()
            visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='hot', new_fig=False, save_fig='psi_flux_tubes%03d' % ctr)
            pl.figure()
            filled_arr = b_mag.copy()
            filled_arr[flux_tube_mask] = 2*b_mag.max()
            visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='hot', new_fig=False, save_fig='b_mag_flux_tubes%03d' % ctr)
            pl.close('all')
            if 0:
                filled_arr = psi_arr.copy()
                def hf((a,b)): return height_func(a), height_func(b)
                ctr = 0
                for region in sorted(regions, key=hf, reverse=True):
                    X = [_[0] for _ in regions[region]]
                    Y = [_[1] for _ in regions[region]]
                    filled_arr[X,Y] = 2*psi_arr.max()
                    if not ctr:
                        visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='gray', new_fig=False)
                        raw_input("enter to continue")
                    ctr += 1
                    ctr %= len(regions) / 20
                visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='gray', new_fig=False)
                raw_input("enter to continue")
                pl.close('all')
        ctr += 1

        del c_tree
        del regions
        del cpts#}}}

from region_analysis import radial_profiles
def test_radial_profiles():
    psi_arr = nth_timeslice('data.h5', 'psi', -1)
    bx_arr = nth_timeslice('data.h5', 'bx', -1)
    by_arr = nth_timeslice('data.h5', 'by', -1)
    # for bx_arr, by_arr, psi_arr in h5_gen('data.h5', ('bx', 'by', 'psi')):
    psi_arr = psi_arr.astype(np.double)
    mask = np.zeros(psi_arr.shape, dtype=np.bool_)
    bmag = np.sqrt(bx_arr**2 + by_arr**2).astype(np.double)
    surf = _cp.TopoSurface(psi_arr)
    rprofs = radial_profiles(surf, threshold=50, expand_regions=10,
            other_arr=bmag, mask=mask)
    # psi_arr[mask] = 2 * psi_arr.max()
    # pl.figure()
    # pl.imshow(psi_arr, interpolation='nearest')
    # pl.figure()
    # pl.imshow(bmag)
    bmag_fig = pl.figure()
    fluxes_vs_dists = pl.figure()
    fluxes_vs_flux = pl.figure()
    pl.hold(False)
    for minmax, (rprof, region) in rprofs.items():
        # minmax_flux = arr_div[minmax]
        rprof.sort(key=lambda x: x[1])
        import pdb; pdb.set_trace()
        pts, fluxes, avg_b, avg_b_errs, avg_dists, avg_dists_errs = zip(*rprof)
        region_xs, region_ys = zip(*region)
        mask = np.ones(bmag.shape, dtype=np.bool_)
        mask[region_xs, region_ys] = False
        masked_bmag = bmag.copy()
        masked_bmag[mask] = 0.0
        pl.figure(bmag_fig.number)
        pl.imshow(masked_bmag)
        avg_b = np.array(avg_b)
        avg_b -= avg_b[0]
        avg_dists = np.array(avg_dists)
        fluxes = np.array(fluxes)
        fluxes -= fluxes.min()
        # fluxes = np.abs(np.array(fluxes) - minmax_flux)
        # avg_fluxes = np.abs(np.array(avg_fluxes) - minmax_flux)
        # pl.plot(avg_dists, avg_fluxes, 'd-')
        b_over_flux = (avg_b-(avg_b.min())) / fluxes
        pl.figure(fluxes_vs_dists.number)
        pl.plot(fluxes, avg_b, 'd-')
        # pl.errorbar(avg_dists, b_over_flux, yerr=avg_fluxes_errs,
                # hold=True)
        pl.figure(fluxes_vs_flux.number)
        # delta_r = np.diff(avg_dists)
        delta_flux = np.diff(fluxes)
        # dr_b_over_r = np.diff(b_over_r) / delta_r
        dr_b_over_flux = np.diff(b_over_flux) / delta_flux
        # pl.plot(avg_dists, (avg_fluxes-avg_fluxes.min()) / avg_dists, 'x-')
        pl.plot(fluxes[:-1], dr_b_over_flux, 'x-')
        # pl.errorbar(avg_dists[:-1], dr_b_over_r, 
                # xerr=avg_dists_errs[:-1],
                # yerr=avg_fluxes_errs[:-1],
                # hold=True)
        raw_input("enter to continue")
    pl.figure(fluxes_vs_dists.number)
    pl.grid()
    pl.figure(fluxes_vs_flux.number)
    pl.grid()
    import pdb; pdb.set_trace()
    pl.close('all')

test_radial_profiles()
