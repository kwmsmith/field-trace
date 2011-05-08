from composite_vis import visualize
from scipy import ndimage
import numpy as np

from field_trace import _critical_points as _cp
from field_trace.upsample import upsample

from nose.tools import eq_, ok_, set_trace

import pylab as pl
pl.ion()

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
    peaks  = crit_pts.peaks
    passes  = crit_pts.passes
    pits  = crit_pts.pits
    non_morse_passes = set()
    disconnected_pps = set()
    for pss in passes:
        cpeaks = snet.predecessors(pss)
        cpits  = snet.successors(pss)
        if len(cpeaks) != 2 or len(cpits) != 2:
            non_morse_passes.add(pss)
        for cpk in cpeaks:
            assert cpk in peaks
        for cpit in cpits:
            assert cpit in pits
    for pp in pits.union(peaks):
        if pp in pits:
            pp_nbrs = snet.predecessors(pp)
        elif pp in peaks:
            pp_nbrs = snet.successors(pp)
        if not pp_nbrs:
            disconnected_pps.add(pp)
        for pn in pp_nbrs:
            assert pn in passes
    return non_morse_passes, disconnected_pps


def test_critical_points():

    def _tester(arr, vis=False):
        surf = _cp.TopoSurface(arr)
        peaks = surf.crit_pts.peaks
        pits = surf.crit_pts.pits
        passes = surf.crit_pts.passes
        print "\npeaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        snet = surf.surf_network
        snet_points = set(snet.nodes())
        missed_passes = passes.difference(snet_points)
        missed_pits = pits.difference(snet_points)
        missed_peaks = peaks.difference(snet_points)
        print "pits + passes = %d" % (len(peaks) + len(pits))
        print "missed passes: %d" % len(missed_passes)
        print "missed peaks: %d" % len(missed_peaks)
        print "missed pits: %d" % len(missed_pits)
        # regions = surf.get_minmax_regions()
        # mask = np.zeros(arr.shape, dtype=np.bool_)
        # for region in regions:
            # for p in region:
                # mask[p] = True
        # mask_arr = arr.copy()
        # mask_arr[mask] = 2 * arr.max()
        # if vis:
            # visualize(mask_arr, mesh=None, crit_pts=surf.crit_pts, surf_network=snet)
            # raw_input('enter to continue')
        non_morse = verify_snet(snet, surf.crit_pts)
        eq_(len(missed_passes), 0)
        eq_(len(missed_peaks), 0)
        eq_(len(missed_pits), 0)
        # eq_(len(peaks) + len(pits), len(passes))
        eq_(len(snet_points), len(peaks)+len(pits)+len(passes))
        surf.simplify_surf_network(measure=surf.get_peak_pit_region_area, threshold=4)
        # regions = surf.get_minmax_regions()
        # mask = np.zeros(arr.shape, dtype=np.bool_)
        # for region in regions:
            # for p in region:
                # mask[p] = True
        # mask_arr = arr.copy()
        # mask_arr[mask] = 2 * arr.max()
        # if vis:
            # visualize(mask_arr, mesh=None, crit_pts=surf.crit_pts, surf_network=snet)
            # raw_input('enter to continue')
        # reeb = surf.get_reeb_graph()

    for _ in range(1):
        yield _tester, random_periodic_upsample(128, 4, seed=_), False


def compare_graphs(g1, g2):
    g1_unordered = dict((k, set(v)) for (k, v) in g1.items())
    g2_unordered = dict((k, set(v)) for (k, v) in g2.items())
    return g1_unordered == g2_unordered

class cpts(object):
    def __init__(self, peaks, pits, passes):
        self.peaks = peaks
        self.pits = pits
        self.passes = passes

def test_reeb1():
    surf_net = _cp.netx.DiGraph()
    node2h_map = {
            6: 0,
            5: 1,
            4: 2,
            3: 3,
            2: 4,
            1: 5,
            }
    node2h = lambda n: (node2h_map[n], n)

    def _add_edge(n1, n2):
        surf_net.add_edge(n1, n2, dh=abs(node2h(n1)[0] - node2h(n2)[0]))

    _add_edge(1,3)
    _add_edge(1,4)
    _add_edge(2,3)
    _add_edge(2,4)
    _add_edge(3,5)
    _add_edge(3,6)
    _add_edge(4,5)
    _add_edge(4,6)
    peaks  = set([1, 2])
    passes = set([3, 4])
    pits   = set([5, 6])
    crit_pts = cpts(peaks=peaks, passes=passes, pits=pits)
    reeb_gr = _cp.get_reeb_graph(surf_net, crit_pts, node2h)
    ok_(
        compare_graphs(
            reeb_gr.edge,
            {1: [3], 2: [3], 3: [1, 2, 4], 4: [3, 5, 6], 5: [4], 6: [4]}))

def test_reeb2():
    surf_net = _cp.netx.DiGraph()
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

    def _add_edge(n1, n2):
        surf_net.add_edge(n1, n2, dh=abs(node2h(n1)[0] - node2h(n2)[0]))

    _add_edge(1,5)
    _add_edge(1,6)
    _add_edge(1,7)
    _add_edge(2,7)
    _add_edge(3,6)
    _add_edge(4,5)
    _add_edge(5,8)
    _add_edge(6,8)
    _add_edge(7,8)
    peaks  = set([1, 2, 3, 4])
    passes = set([5, 6, 7])
    pits   = set([8])
    crit_pts = cpts(peaks=peaks, passes=passes, pits=pits)
    reeb_gr = _cp.get_reeb_graph(surf_net, crit_pts, node2h)
    ok_(
        compare_graphs(
            reeb_gr.edge,
            {1: [7], 2: [7], 3: [6], 4: [5], 5: [4, 8, 6], 6: [3, 7, 5], 7: [1, 2, 6], 8: [5]}))
    return

def _test_reeb3():
    surf_net = _cp.netx.DiGraph()
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

