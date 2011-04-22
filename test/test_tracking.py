from os import path
import tables
import tracking
from tracking import TrackRegion
import numpy as np
import _critical_points as _cp
from cPickle import load, dump

from contour_tree import wraparound_dist

import pylab as pl

from nose.tools import ok_, eq_
from scipy.ndimage.filters import gaussian_filter

def get_cached_surfs(h5fname, sigma=None, force=False):
    saved_name = path.realpath(h5fname)+'_pickled_snets'
    if not force:
        try:
            fh = open(saved_name, 'rb')
        except IOError:
            pass
        else:
            surfs = load(fh)
            sigma = load(fh)
            fh.close()
            return surfs, sigma
    # rebuild the regions list.
    dta = tables.openFile(h5fname, 'r')
    psigp = dta.getNode('/psi')._v_children
    surfs = []
    for psi_arr_name in sorted(psigp):
        print psi_arr_name
        psi_arr = psigp[psi_arr_name].read()
        if sigma is not None:
            psi_arr = gaussian_filter(psi_arr, sigma=sigma, mode='wrap')
        surf = _cp.TopoSurface(psi_arr)
        surfs.append(surf)
    dta.close()
    with open(saved_name, 'w') as fh:
        dump(surfs, fh)
        dump(sigma, fh)
    return surfs, sigma

def extract_regions(h5fname, sigma=None):
    dta = tables.openFile(h5fname, 'r')
    psigp = dta.getNode('/psi')._v_children
    regions = []
    for idx, psi_arr_name in enumerate(sorted(psigp)):
        print psi_arr_name
        psi_arr = psigp[psi_arr_name].read()
        if sigma is not None:
            psi_arr = gaussian_filter(psi_arr, sigma=sigma, mode='wrap')
        surf = _cp.TopoSurface(psi_arr)
        # surf.simplify_surf_network(surf.get_peak_pit_region_area, threshold=100)
        psi_regions = surf.get_minmax_regions()
        tslice = []
        for (cpt, type) in psi_regions:
            area = len(psi_regions[cpt, type])
            minmax = psi_arr[cpt]
            tslice.append(TrackRegion(loc=cpt, area=area, val=minmax, ispeak=(type=='peak')))
        regions.append((idx, tslice))
    dta.close()
    return regions

def extract_regions_cached(h5fname, sigma=None, force=False):
    saved_name = path.realpath(h5fname)+'_pickled_regions'
    if not force:
        try:
            fh = open(saved_name, 'rb')
        except IOError:
            pass
        else:
            regions = load(fh)
            sigma = load(fh)
            fh.close()
            return regions, sigma
    # rebuild the regions list.
    regions = extract_regions(h5fname, sigma=sigma)
    with open(saved_name, 'w') as fh:
        dump(regions, fh)
        dump(sigma, fh)
    return regions, sigma

def track_regions_image(slice_pts, h5fname, sigma):
    tslice_to_regions = extract_regions('data.h5', sigma=sigma)
    dta = tables.openFile(h5fname, 'r')
    psigp = dta.getNode('/psi')._v_children
    pl.ioff()
    for idx, psi_arr_name in enumerate(sorted(psigp)):
        print idx
        psi_arr = psigp[psi_arr_name].read()
        pl.figure()
        pl.imshow(psi_arr, interpolation='nearest', cmap='hot')
        for (loc, num) in slice_pts[idx]:
            pl.text(loc[1], loc[0], str(num))
        xs = [reg.loc[0] for reg in tslice_to_regions[idx][1]]
        ys = [reg.loc[1] for reg in tslice_to_regions[idx][1]]
        pl.scatter(ys, xs, c='r', marker='s')
        pl.savefig('psi_%03d.pdf' % idx)
        pl.close('all')
    dta.close()

def plot_track_props(tracks, nx, ny):
    pl.ioff()
    wdist = wraparound_dist(nx, ny)
    val_fig = pl.figure()
    area_fig = pl.figure()
    psn_fig = pl.figure()
    delta_vals = []
    delta_dists = []
    for tr in tracks:
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


def test_track_greedy():
    nx, ny = 512, 512
    sigma = 4.0
    tslice_to_regions = extract_regions('data.h5', sigma=sigma)
    tracks = tracking.track_regions_greedy(tslice_to_regions, nx, ny)
    reverse_tracks = tracking.track_regions_greedy(list(reversed(tslice_to_regions)), nx, ny, reverse=True)
    pl.ion()
    # val_fig = pl.figure()
    # area_fig = pl.figure()
    # psn_fig = pl.figure()
    maxlen = len(tslice_to_regions)
    # tracks = [tr for tr in tracks if len(tr) >= 0.3 * maxlen]
    # rev_tracks = [tr for tr in reverse_tracks if len(tr) >= 0.3 * maxlen]
    tracks = reverse_tracks
    from collections import defaultdict
    slice_pts = defaultdict(list)
    for num, tr in enumerate(tracks):
        for idx, reg in tr:
            slice_pts[idx].append((reg.loc, num))
    plot_track_props(tracks, nx, ny)
    track_regions_image(slice_pts, 'data.h5', sigma=sigma)
    # for tr in tracks:
        # idxs, regs = zip(*tr)
        # pl.figure(val_fig.number)
        # pl.plot(idxs, [reg.val for reg in regs], 's-', hold=True)
        # pl.figure(area_fig.number)
        # pl.semilogy(idxs, [reg.area for reg in regs], 's-', hold=True)
        # pl.figure(psn_fig.number)
        # dists = [wdist(regs[i].loc, regs[i+1].loc) for i in range(len(regs)-1)]
        # pl.plot(idxs[:-1], dists, 's-', hold=True)
    # import pdb; pdb.set_trace()

class _test_track(object):

    def setUp(self):
        self.npts = 500
        self.nx = self.ny = 100
        self.pts1 = np.random.randint(self.nx, size=(self.npts, 3))
        self.pts1[:,2] = np.random.randint(1, 10, size=self.npts)
        self.pts2 = np.random.randint(self.nx, size=(self.npts, 3))
        self.pts2[:,2] = np.random.randint(1, 10, size=self.npts)

    def _test_track_forward_backward(self):
        tracking.track_forward_backward(self.pts1, self.pts2, self.nx, self.ny)


    def test_track(self):
        map = tracking.track_regions(self.pts1, self.pts2, self.nx, self.ny)
        eq_(map, range(self.npts))
        self.pts2 = self.pts1.copy()
        self.pts2[0], self.pts2[1] = self.pts1[1], self.pts1[0]
        map = tracking.track_regions(self.pts1, self.pts2, self.nx, self.ny)
        check = range(self.npts)
        check[0], check[1] = check[1], check[0]
        eq_(map, check)
