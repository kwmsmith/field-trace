import scipy as sp
from scipy import ndimage
import numpy as np

import _critical_points as _cp
from upsample import upsample

from kaw_analysis import test_vcalc as tv

from nose.tools import eq_, ok_, set_trace

def random_periodic_upsample(N, upsamp, seed=None):
    if seed is not None:
        np.random.seed(seed)
    arr = np.random.normal(size=(N/upsamp, N/upsamp))
    return upsample(arr, factor=upsamp)

def random_periodic(N, sigma, seed=None):
    if seed is not None:
        np.random.seed(seed)
    arr = np.random.normal(size=(N,N))
    arr -= arr.mean()
    return ndimage.gaussian_filter(arr, sigma=sigma, mode='wrap')

def add_noise(arr, rms_frac=0.001):
    rms = arr.std()
    noise = np.random.normal(scale=rms_frac*rms, size=arr.shape)
    return arr + noise

def test_critical_points():

    def _tester(arr, vis=False):
        surf = _cp.TopoSurface(arr)
        peaks = surf.crit_pts['peaks']
        pits = surf.crit_pts['pits']
        passes = surf.crit_pts['passes']
        print "\npeaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        snet = surf.surf_network
        snet_points = set(snet._g.keys())
        missed_passes = passes.difference(snet_points)
        missed_pits = pits.difference(snet_points)
        missed_peaks = peaks.difference(snet_points)
        print "missed passes: %d" % len(missed_passes)
        print "missed peaks: %d" % len(missed_peaks)
        print "missed pits: %d" % len(missed_pits)
        eq_(len(missed_passes), 0)
        eq_(len(missed_peaks), 0)
        eq_(len(missed_pits), 0)
        eq_(len(peaks) + len(pits), len(passes))
        eq_(len(snet_points), len(peaks)+len(pits)+len(passes))
        if vis:
            visualize(arr, mesh=None, crit_pts=surf.crit_pts, surf_network=snet)
            raw_input('enter to continue')

    for _ in range(1):
        yield _tester, random_periodic_upsample(64, 8, seed=None), True

def visualize(arr, mesh=None, crit_pts=None, surf_network=None):
    import pylab as pl
    pl.ioff()
    fig = pl.figure()
    # pl.imshow(arr, interpolation='nearest', cmap='jet')
    pl.imshow(arr, cmap='jet', interpolation='nearest')
    nx, ny = arr.shape
    if mesh is not None:
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                other_pts = mesh._g[i,j]
                for x,y in other_pts:
                    pl.plot([j,y], [i,x], 'k--')
    if surf_network is not None:
        for node in surf_network._g:
            node_x, node_y = node
            nbrs = surf_network._g[node]
            for nbr in nbrs:
                nbr_x, nbr_y = nbr
                pl.plot([node_y, nbr_y], [node_x, nbr_x], 'k--')
    if crit_pts is not None:
        pits, passes, peaks = crit_pts['pits'], crit_pts['passes'], crit_pts['peaks']
        X = [_[0] for _ in pits]
        Y = [_[1] for _ in pits]
        pl.scatter(Y, X, marker='o', c='b', s=50)
        X = [_[0] for _ in peaks]
        Y = [_[1] for _ in peaks]
        pl.scatter(Y, X, marker='d', c='r', s=50)
        X = [_[0] for _ in passes]
        Y = [_[1] for _ in passes]
        pl.scatter(Y, X, marker='d', c='k', s=50)
    pl.ion()
    pl.figure(fig.number)
