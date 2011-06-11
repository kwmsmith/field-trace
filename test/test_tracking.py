import numpy as np
import tables
import tracking

from tracking import extract_regions

from contour_tree import wraparound_dist

import pylab as pl

# h5fname = '/Users/ksmith/Research/thesis/data/high-time-res/data.h5'
h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10/data.h5'

def track_regions_image(slice_regions, h5fname, sigma):
    tslice_to_regions = extract_regions(h5fname, sigma=sigma)
    dta = tables.openFile(h5fname, 'r')
    pl.ioff()
    for idx, psi_arr in enumerate(dta.walkNodes('/psi', 'Array')):
        print idx
        psi_arr = psi_arr.read()
        psi_arr_mask = psi_arr.copy()
        for (reg, num) in slice_regions[idx]:
            psi_arr_mask[zip(*reg.region)] = 2*psi_arr.max()
        pl.figure()
        pl.imshow(psi_arr_mask, interpolation='nearest', cmap='hot')
        for (reg, num) in slice_regions[idx]:
            loc = reg.loc
            pl.text(loc[1], loc[0], str(num))
        xs = [reg.loc[0] for reg in tslice_to_regions[idx][1]]
        ys = [reg.loc[1] for reg in tslice_to_regions[idx][1]]
        pl.scatter(ys, xs, c='r', marker='s')
        pl.savefig('psi_%03d.pdf' % idx)
        pl.close('all')
    dta.close()

def plot_track_props(tracks, nx, ny, len_cutoff=20):
    pl.ioff()
    wdist = wraparound_dist(nx, ny)
    val_fig = pl.figure()
    area_fig = pl.figure()
    psn_fig = pl.figure()
    delta_vals = []
    delta_dists = []
    for tr in tracks:
        if len(tr) < len_cutoff: continue
        idxs, regs = zip(*tr)
        delta_vals.extend([abs(regs[idx].val-regs[idx+1].val) for idx in range(len(regs)-1)])
        dists = [wdist(regs[i].loc, regs[i+1].loc) for i in range(len(regs)-1)]
        delta_dists.extend([abs(dists[idx]-dists[idx+1]) for idx in range(len(dists)-1)])
        pl.figure(val_fig.number)
        pl.plot(idxs, [reg.val for reg in regs], 's-', hold=True)
        pl.figure(area_fig.number)
        pl.semilogy(idxs, [reg.area for reg in regs], 's-', hold=True)
        pl.figure(psn_fig.number)
        pl.plot(idxs[:-1], dists, 's-', hold=True)
    pl.figure(val_fig.number)
    pl.savefig("val_v_time.pdf")
    pl.figure(area_fig.number)
    pl.savefig("area_v_time.pdf")
    pl.figure(psn_fig.number)
    pl.savefig("psn_v_time.pdf")
    pl.figure()
    pl.hist(delta_vals, bins=pl.sqrt(len(delta_vals)))
    pl.savefig("delta_vals.pdf")
    pl.figure()
    pl.hist(delta_dists, bins=pl.sqrt(len(delta_dists)))
    pl.savefig("delta_dists.pdf")
    pl.close('all')

def get_slice_regions(tracks):
    from collections import defaultdict
    slice_regions = defaultdict(list)
    tracks = [tr for tr in tracks if len(tr) > 2]
    for num, tr in enumerate(tracks):
        for idx, reg in tr:
            slice_regions[idx].append((reg, num))
    return slice_regions

def test_track_greedy():
    nx, ny = 512, 512
    sigma = None
    tslice_to_regions = extract_regions(h5fname, sigma=sigma)

    tracks = tracking.track_regions_greedy(tslice_to_regions, nx, ny)
    tracks = tracking.chop_tracks(tracks, area_frac=0.1)
    slice_regions = get_slice_regions(tracks)
    track_regions_image(slice_regions, h5fname, sigma=sigma)
