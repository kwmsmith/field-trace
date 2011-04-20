from os import path
import tables
import tracking
from tracking import TrackRegion
import numpy as np
import _critical_points as _cp
from pickle import load, dump

from contour_tree import wraparound_dist

import pylab as pl

from nose.tools import ok_, eq_

def surf_networks(h5fname, force=False):
    saved_name = path.realpath(h5fname)+'_pickled_surf_networks'
    if not force:
        try:
            fh = open(saved_name, 'rb')
        except IOError:
            pass
        else:
            snets = load(fh)
            fh.close()
            return snets
    # rebuild the surface networks
    dta = tables.openFile(h5fname, 'r')
    psigp = dta.getNode('/psi')._v_children
    snets = []
    for psi_arr_name in sorted(psigp)[:1]:
        psi_arr = psigp[psi_arr_name].read()
        surf = _cp.TopoSurface(psi_arr)
        snets.append(surf)
    with open(saved_name, 'w') as fh:
        dump(snets, fh)
    return snets

def extract_regions(h5fname, force=False):
    snets = surf_networks(h5fname, force=force)
    regions = []
    for snet in snets:
        psi_regions = snet.get_minmax_regions()
        tslice = []
        for (cpt, type) in psi_regions:
            area = len(psi_regions[cpt, type])
            minmax = snet.arr[cpt]
            tslice.append(TrackRegion(loc=cpt, area=area, val=minmax, ispeak=(type=='peak')))
        regions.append(tslice)
    return regions

def track_regions_image(slice_pts, h5fname):
    all_regions = extract_regions('data.h5')
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
        xs = [reg.loc[0] for reg in all_regions[idx]]
        ys = [reg.loc[1] for reg in all_regions[idx]]
        pl.scatter(ys, xs, c='r', marker='s')
        pl.savefig('psi_%03d.pdf' % idx)
        pl.close('all')
    dta.close()

def test_track_greedy():
    nx, ny = 512, 512
    wdist = wraparound_dist(nx, ny)
    all_regions = extract_regions('data.h5')
    tracks = tracking.track_regions_greedy(all_regions, nx, ny)
    pl.ion()
    # val_fig = pl.figure()
    # area_fig = pl.figure()
    # psn_fig = pl.figure()
    maxlen = len(all_regions)
    tracks = [tr for tr in tracks if len(tr) >= 0.5 * maxlen]
    from collections import defaultdict
    slice_pts = defaultdict(list)
    for num, tr in enumerate(tracks):
        for idx, reg in tr:
            slice_pts[idx].append((reg.loc, num))
    track_regions_image(slice_pts, 'data.h5')
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
