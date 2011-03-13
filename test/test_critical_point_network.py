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
        print "\nmeshing array"
        gr = _cp.mesh(arr)
        print "classifying nodes"
        classes = _cp.classify_nodes(arr, gr)
        passes = classes['passes']
        pits = classes['pits']
        peaks = classes['peaks']
        eq_(len(peaks) + len(pits), len(passes))
        print "peaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        print "getting surface network"
        snet = _cp.surface_network(arr, gr, passes, peaks, pits)
        snet_points = set(snet._g.keys())
        missed_passes = passes.difference(snet_points)
        missed_pits = pits.difference(snet_points)
        missed_peaks = peaks.difference(snet_points)
        print "missed passes: %d" % len(missed_passes)
        print "missed peaks: %d" % len(missed_peaks)
        print "missed pits: %d" % len(missed_pits)
        if vis:
            visualize(arr, gr=None, classes=classes, surf_network=None)
            raw_input('enter to continue')

    for _ in range(4):
        yield _tester, random_periodic_upsample(128, 4, seed=None), True

def visualize(arr, gr=None, classes=None, surf_network=None):
    import pylab as pl
    pl.ioff()
    fig = pl.figure()
    # pl.imshow(arr, interpolation='nearest', cmap='jet')
    pl.imshow(arr, cmap='jet', interpolation='nearest')
    nx, ny = arr.shape
    if gr is not None:
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                other_pts = gr._g[i,j]
                for x,y in other_pts:
                    pl.plot([j,y], [i,x], 'k--')
    if surf_network is not None:
        for node in surf_network._g:
            node_x, node_y = node
            nbrs = surf_network._g[node]
            for nbr in nbrs:
                nbr_x, nbr_y = nbr
                pl.plot([node_y, nbr_y], [node_x, nbr_x], 'k--')
    if classes is not None:
        pits, passes, peaks = classes['pits'], classes['passes'], classes['peaks']
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
