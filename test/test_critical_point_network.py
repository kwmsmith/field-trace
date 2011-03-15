import scipy as sp
from scipy import ndimage
import numpy as np

import _critical_points as _cp
from upsample import upsample

from kaw_analysis import test_vcalc as tv

from nose.tools import eq_, ok_, set_trace

from pprint import pprint

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

def verify_snet(snet, crit_pts):
    peaks  = crit_pts['peaks']
    passes  = crit_pts['passes']
    pits  = crit_pts['pits']
    for pss in passes:
        pass_nbrs = snet._g[pss]
        eq_(len(pass_nbrs), 4)
        for pn in pass_nbrs:
            ok_(pn in peaks or pn in pits)
    for pp in pits.union(peaks):
        pp_nbrs = snet._g[pp]
        for pn in pp_nbrs:
            ok_(pn in passes)


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
        print "pits + passes = %d" % (len(peaks) + len(pits))
        print "missed passes: %d" % len(missed_passes)
        print "missed peaks: %d" % len(missed_peaks)
        print "missed pits: %d" % len(missed_pits)
        eq_(len(missed_passes), 0)
        eq_(len(missed_peaks), 0)
        eq_(len(missed_pits), 0)
        eq_(len(peaks) + len(pits), len(passes))
        eq_(len(snet_points), len(peaks)+len(pits)+len(passes))
        # verify_snet(snet, surf.crit_pts)
        reeb = surf.get_reeb_graph()
        if vis:
            visualize(arr, mesh=None, crit_pts=surf.crit_pts, surf_network=snet)
            raw_input('enter to continue')

    for _ in range(10):
        yield _tester, random_periodic_upsample(32, 8, seed=_), False

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

def compare_graphs(g1, g2):
    g1_unordered = dict((k, set(v)) for (k, v) in g1.items())
    g2_unordered = dict((k, set(v)) for (k, v) in g2.items())
    return g1_unordered == g2_unordered

def test_reeb1():
    surf_net = _cp.graph()
    node2h_map = {
            6: 0,
            5: 1,
            4: 2,
            3: 3,
            2: 4,
            1: 5,
            }
    node2h = lambda n: (node2h_map[n], n)
    surf_net.add_edge(1,3)
    surf_net.add_edge(1,4)
    surf_net.add_edge(2,3)
    surf_net.add_edge(2,4)
    surf_net.add_edge(3,5)
    surf_net.add_edge(3,6)
    surf_net.add_edge(4,5)
    surf_net.add_edge(4,6)
    peaks  = set([1, 2])
    passes = set([3, 4])
    pits   = set([5, 6])
    crit_pts = {'peaks': peaks,
                'passes': passes,
                'pits' : pits}
    reeb_gr = _cp.get_reeb_graph(surf_net, crit_pts, node2h)
    ok_(
        compare_graphs(
            reeb_gr._g,
            {1: [3], 2: [3], 3: [1, 2, 4], 4: [3, 5, 6], 5: [4], 6: [4]}))

def test_reeb2():
    surf_net = _cp.graph()
    node2h_map = {
            8: 0,
            7: 1,
            6: 1,
            5: 1,
            4: 2,
            3: 2,
            2: 2,
            1: 3,
            }
    node2h = lambda n: (node2h_map[n], n)
    surf_net.add_edge(1,5)
    surf_net.add_edge(1,6)
    surf_net.add_edge(1,7)
    surf_net.add_edge(2,7)
    surf_net.add_edge(3,6)
    surf_net.add_edge(4,5)
    surf_net.add_edge(5,8)
    surf_net.add_edge(6,8)
    surf_net.add_edge(7,8)
    peaks  = set([1, 2, 3, 4])
    passes = set([5, 6, 7])
    pits   = set([8])
    crit_pts = {'peaks': peaks,
                'passes': passes,
                'pits' : pits}
    reeb_gr = _cp.get_reeb_graph(surf_net, crit_pts, node2h)
    ok_(
        compare_graphs(
            reeb_gr._g,
            {1: [7], 2: [7], 3: [6], 4: [5], 5: [4, 8, 6], 6: [3, 7, 5], 7: [1, 2, 6], 8: [5]}))
    return

def test_reeb3():
    surf_net = _cp.graph()
    surf_net._g = {
            (5, 15): [(22, 28), (6, 15)],
            (5, 31): [(6, 25), (26, 6), (10, 2)],
            (6, 15): [(27, 0), (5, 15), (6, 16)],
            (6, 16): [(6, 25), (14, 9), (6, 15)],
            (6, 25): [(27, 0), (21, 18), (6, 16), (5, 31)],
            (10, 2): [(21, 18), (27, 0), (5, 31), (16, 3)],
            (14, 9): [(27, 0), (21, 18), (16, 3), (6, 16)],
            (16, 3): [(22, 28), (14, 9), (26, 6), (10, 2)],
            (21, 18): [(6, 25), (22, 28), (14, 9), (26, 6), (10, 2)],
            (22, 28): [(21, 18), (27, 0), (5, 15), (16, 3)],
            (26, 6): [(27, 0), (21, 18), (16, 3), (5, 31)],
            (27, 0): [(6, 25), (22, 28), (14, 9), (26, 6), (6, 15), (10, 2)],
            }
    crit_pts =  {
            'passes': set([(6, 15), (6, 25), (10, 2), (14, 9), (22, 28), (26, 6)]),
            'peaks': set([(21, 18), (27, 0)]),
            'pits': set([(5, 15), (5, 31), (6, 16), (16, 3)]),
            }
    node2h_map = {
            (5, 31): (-0.014356736424827819, (5, 31)),
            (5, 15): (-0.027839962424853988, (5, 15)),
            (6, 25): (-0.0055548455026826482, (6, 25)),
            (6, 15): (-0.027593845591386536, (6, 15)),
            (6, 16): (-0.027642410602253194, (6, 16)),
            (10, 2): (-0.0075859636696798422, (10, 2)),
            (14, 9): (-0.0072652858280034555, (14, 9)),
            (16, 3): (-0.012332851686022498, (16, 3)),
            (21, 18): (0.039444128294746132, (21, 18)),
            (22, 28): (0.0016904116601130521, (22, 28)),
            (26, 6): (0.0064175915669918903, (26, 6)),
            (27, 0): (0.011876681723300701, (27, 0)),
            }
    node2h = lambda n: node2h_map[n]
    reeb_gr = _cp.get_reeb_graph(surf_net, crit_pts, node2h)
    set_trace()
