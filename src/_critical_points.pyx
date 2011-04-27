# cython: profile=True

cimport numpy as np
import numpy as np
from pprint import pprint, pformat

import networkx as netx

import heapq

cdef int AC = 0
cdef int BD = 1

cdef int envelope(double a, double b, double c, double d):
    if b + d >= a + c:
        return BD
    else:
        return AC

cdef int flattest(double a, double b, double c, double d):
    if (d - b)*(d - b) >= (a - c)*(a - c):
        return AC
    else:
        return BD

cdef int connect_diagonal(double a, double b, double c, double d):
    """
    returns 'ac' or 'bd' indicating which nodes to connect in the square cell.
    
    """
    return flattest(a, b, c, d)
    # return envelope(a, b, c, d)
    # return AC

DEF NNBRS = 6
cdef void _neighbors6(int x, int y, int nx, int ny, int *result):
    '''
    result is a 12 element integer array that holds the neighbor's x & y components.
    '''
    result[0] = (x-1+nx) % nx
    result[1] = (y-1+ny) % ny

    result[2] = (x-1+nx) % nx
    result[3] = y

    result[4] = x
    result[5] = (y+1) % ny

    result[6] = (x+1) % nx
    result[7] = (y+1) % ny

    result[8] = (x+1) % nx
    result[9] = y

    result[10] = x
    result[11] = (y-1+ny) % ny

def neighbors6(int x, int y, int nx, int ny):
    cdef np.ndarray[np.int32_t, ndim=2] result = np.empty((NNBRS,2), dtype=np.int32)
    _neighbors6(x, y, nx, ny, <int*>result.data)
    return [tuple(r) for r in result]

class graph(object):

    def __init__(self):
        from collections import defaultdict
        self._g = defaultdict(list)

    def add_edge(self, a, b):
        if b not in self._g[a]:
            self._g[a].append(b)
            self._g[b].append(a)

    def remove_node(self, nd):
        if nd in self._g:
            del self._g[nd]

    def remove_edge(self, a, b):
        if b in self._g[a]:
            self._g[a].remove(b)
            self._g[b].remove(a)

    def order_neighbors(self, shape):
        nx, ny = shape
        for node in self._g:
            i, j = node
            neighbors = self._g[node]
            ord_neighbors = []
            for n in (((i-1)%nx, (j-1)%ny),
                      (i,        (j-1)%ny),
                      ((i+1)%nx, (j-1)%ny),
                      ((i+1)%nx,  j      ),
                      ((i+1)%nx, (j+1)%ny),
                      (i       , (j+1)%ny),
                      ((i-1)%nx, (j+1)%ny),
                      ((i-1)%nx,  j      )):
                if n in neighbors:
                    ord_neighbors.append(n)
            assert 4 <= len(neighbors) <= 8
            assert set(ord_neighbors) == set(neighbors)
            self._g[node] = ord_neighbors

    def order_by_key(self, keyfunc):
        for node in self._g:
            self._g[node].sort(key=keyfunc)

    def deepcopy(self):
        from copy import deepcopy
        g = graph()
        g._g = deepcopy(self._g)
        return g

def sort_by_h(gr, arr):
    def _keyfunc(n):
        return (arr[n], n)
    sgr = {}
    for node in gr._g:
        h_sort = sorted(gr._g[node], key=_keyfunc)
        sgr[node] = h_sort
    return sgr

cdef _get_highest_neighbor(np.ndarray[double, ndim=2] arr, node):
    cdef int nbrs[2*NNBRS]
    cdef int i, highest_x, highest_y
    nx = arr.shape[0]
    ny = arr.shape[1]
    _neighbors6(node[0], node[1], nx, ny, nbrs)
    cdef double highest_val = -1e100, val
    for i in range(NNBRS):
        val = arr[nbrs[2*i], nbrs[2*i+1]]
        if val > highest_val:
            highest_val = val
            highest_x = nbrs[2*i]
            highest_y = nbrs[2*i+1]
    return highest_x, highest_y

cdef _get_lowest_neighbor(np.ndarray[double, ndim=2] arr, node):
    cdef int nbrs[2*NNBRS]
    cdef int i, lowest_x, lowest_y
    nx = arr.shape[0]
    ny = arr.shape[1]
    _neighbors6(node[0], node[1], nx, ny, nbrs)
    cdef double lowest_val = 1e100, val
    for i in range(NNBRS):
        val = arr[nbrs[2*i], nbrs[2*i+1]]
        if val < lowest_val:
            lowest_val = val
            lowest_x = nbrs[2*i]
            lowest_y = nbrs[2*i+1]
    return lowest_x, lowest_y

from collections import namedtuple
crit_pts = namedtuple('crit_pts', 'peaks pits passes')

class TopoSurface(object):

    def __init__(self, arr):
        self.nx, self.ny = arr.shape
        self.arr = np.asarray(arr, dtype=np.double)
        self.crit_pts = self.get_crit_pts()
        self._surf_network = None

    # def __getstate__(self):
        # return {
                # 'arr': self.arr,
                # 'mesh': self.mesh,
                # 'h_sorted_mesh': self.h_sorted_mesh,
                # 'peaks': self.crit_pts.peaks,
                # 'pits' : self.crit_pts.pits,
                # 'passes': self.crit_pts.passes,
                # }

    # def __setstate__(self, state):
        # self.arr = state['arr']
        # self.mesh = state['mesh']
        # self.h_sorted_mesh = state['h_sorted_mesh']
        # self.crit_pts = crit_pts(peaks=state['peaks'],
                                 # pits=state['pits'],
                                 # passes=state['passes'])
        # self._surf_network = None

    def get_minmax_regions(self):
        passes = self.crit_pts.passes
        pit_regions = dict([((pit, 'pit'), self.get_minmax_region(pit, passes)) \
                for pit in self.crit_pts.pits])
        peak_regions = dict([((peak, 'peak'), self.get_minmax_region(peak, passes)) \
                for peak in self.crit_pts.peaks])
        all_regions = pit_regions
        all_regions.update(peak_regions)
        return all_regions

    def get_minmax_region(self, node, passes):
        if node in self.crit_pts.pits:
            sign = 1
        elif node in self.crit_pts.peaks:
            sign = -1
        else:
            raise ValueError("node not in peaks or pits")
        cdef set region = set([node])
        cdef list frontier = []
        cdef int nbrs[NNBRS*2]
        _neighbors6(node[0], node[1], self.nx, self.ny, nbrs)
        for i in range(NNBRS):
            n = nbrs[i*2], nbrs[(2*i)+1]
            frontier.append((self.node_height(n, sign=sign), n))
        cdef set frontier_set = set([n for h,n in frontier])
        heapq.heapify(frontier)
        while True:
            h, n = heapq.heappop(frontier)
            frontier_set.remove(n)
            region.add(n)
            if n in passes:
                break
            _neighbors6(n[0], n[1], self.nx, self.ny, nbrs)
            for i in range(NNBRS):
                nbr = nbrs[2*i], nbrs[2*i+1]
                if nbr not in region and nbr not in frontier_set:
                    frontier_set.add(nbr)
                    heapq.heappush(frontier, (self.node_height(nbr, sign=sign), nbr))
        return region

    def get_crit_pts(self):
        cdef np.ndarray[double, ndim=2] arr = self.arr
        cdef int nbrs[2*NNBRS]
        cdef double diffs[NNBRS], node_val
        cdef double diff_neg, diff_pos, diff1, diff2
        cdef int n_change, idx, ix, iy, ii
        cdef set peaks = set()
        cdef set pits = set()
        cdef set passes = set()
        for ix in range(self.nx):
            for iy in range(self.ny):
                _neighbors6(ix, iy, self.nx, self.ny, nbrs)
                node_val = arr[ix, iy]
                node = (ix, iy)
                for ii in range(NNBRS):
                    diffs[ii] = arr[nbrs[2*ii], nbrs[2*ii+1]] - node_val
                diff_neg = diff_pos = 0.0
                for ii in range(NNBRS):
                    dd = diffs[ii]
                    if dd < 0:
                        diff_neg += -dd
                    else:
                        diff_pos += dd
                n_change = 0
                for idx in range(NNBRS):
                    diff1 = diffs[(idx-1+NNBRS)%NNBRS]
                    diff2 = diffs[idx]
                    # if diffs[idx] == 0.0 and diffs[idx-1] * diffs[(idx+1)%len(diffs)] < 0:
                        # n_change += 1
                    # elif diffs[idx-1] * diffs[idx] < 0:
                    if diff1 * diff2 < 0:
                        n_change += 1
                if n_change == 0:
                    if diff_neg > 0 and diff_pos == 0.0:
                        peaks.add(node)
                    elif diff_pos > 0 and diff_neg == 0.0:
                        pits.add(node)
                elif (n_change >= 4 and not n_change % 2) and (diff_pos + diff_neg) > 0:
                    passes.add(node)
        return crit_pts(peaks=peaks, pits=pits, passes=passes)
    
    def nearest_extrema(self, pss):
        cdef int pass_nbrs[2*NNBRS]
        _neighbors6(pss[0], pss[1], self.nx, self.ny, pass_nbrs)
        peaks_n_pits = self.crit_pts.peaks.union(self.crit_pts.pits)
        connected_peaks = set()
        connected_pits = set()
        for i in range(NNBRS):
            nbr = pass_nbrs[2*i], pass_nbrs[2*i+1]
            # go up from every neighbor.
            cur = nbr
            while True:
                if cur in self.crit_pts.peaks:
                    connected_peaks.add(cur)
                    break
                cur = _get_highest_neighbor(self.arr, cur)
                if cur == pss or cur in pass_nbrs:
                    # returned to a previously visited spot; ignore these.
                    break
            # go down from every neighbor.
            cur = nbr
            while True:
                if cur in self.crit_pts.pits:
                    connected_pits.add(cur)
                    break
                cur = _get_lowest_neighbor(self.arr, cur)
                if cur == pss or cur in pass_nbrs:
                    # returned to a previously visited spot; ignore these.
                    break
        return (connected_peaks, connected_pits)

    def simplify_surf_network(self, measure, threshold):
        import time
        snet = self.surf_network
        self.weight_peaks_n_pits(measure)
        pqueue = [(m, node) for (node, m) in snet.nodes(data=True) if m != {}]
        heapq.heapify(pqueue)
        seen = set()
        while True:
            m, node = heapq.heappop(pqueue)
            if node in seen:
                continue
            seen.add(node)
            if m > threshold:
                break
            changed_nodes = self.contract_surf_network(node, measure=measure)
            for cn in changed_nodes:
                heapq.heappush(pqueue, (snet.node[cn], cn))

    def is_removable_edge(self, snet, pss, node, peakorpit):
        if peakorpit == 'peak':
            nbrs = snet.successors(node)
            if len(nbrs) == 1:
                return False
            nearest = max([self.node_height(p) for p in nbrs])
        elif peakorpit == 'pit':
            nbrs = snet.predecessors(node)
            if len(nbrs) == 1:
                return False
            nearest = min([self.node_height(p) for p in nbrs])
        if nearest[1] == pss:
            return False
        return True

    def is_morse(self, pss):
        snet = self.surf_network
        if len(snet.predecessors(pss)) != 2:
            return False
        if len(snet.successors(pss)) != 2:
            return False
        return True

    def is_orphan_peak_pit(self, snet, pp):
        return not (snet.successors(pp) + snet.predecessors(pp))

    def regularize_surf_network(self, snet):
        crit_pts = self.crit_pts
        # every pass should be connected to 2 distinct peaks & 2 distinct pits.
        for pss in crit_pts.passes:
            cpeaks = snet.predecessors(pss)
            cpits  = snet.successors(pss)
            if len(cpeaks) > 2:
                onodes = sorted([(self.node_height(cpeak), cpeak) for cpeak in cpeaks])
                for h, onode in onodes:
                    if len(snet.predecessors(pss)) == 2:
                        break
                    if self.is_removable_edge(snet, pss, onode, 'peak'):
                        snet.remove_edge(onode, pss)
                        assert not self.is_orphan_peak_pit(snet, onode)
            if len(cpits) > 2:
                onodes = sorted([(self.node_height(cpit), cpit) for cpit in cpits])
                onodes.reverse()
                for h, onode in onodes:
                    if len(snet.successors(pss)) == 2:
                        break
                    if self.is_removable_edge(snet, pss, onode, 'pit'):
                        snet.remove_edge(pss, onode)
                        assert not self.is_orphan_peak_pit(snet, onode)
            assert not self.is_orphan_peak_pit(snet, pss)
            # assert len(snet.predecessors(pss)) <= 2
            # assert len(snet.successors(pss)) <= 2

    def contract_surf_network(self, node, measure=None):
        snet = self.surf_network
        if node in self.crit_pts.pits:
            passes = snet.pred[node]
        elif node in self.crit_pts.peaks:
            passes = snet.succ[node]
        else:
            raise ValueError("given node not in peaks or pits")
        passes = [(dd['dh'], p) for p, dd in passes.items()]
        passes.sort()
        if not passes:
            return set()
        c_pass_h, c_pass = passes[0]
        if not self.is_morse(c_pass):
            return set()
        rest_passes = passes[1:]
        if node in self.crit_pts.pits:
            other_nodes = snet.succ[c_pass]
        elif node in self.crit_pts.peaks:
            other_nodes = snet.pred[c_pass]
        other_nodes = [(dd['dh'], p) for p, dd in other_nodes.items() if p != node]
        other_nodes.sort()
        to_update = set()
        for pss_dh, pss in rest_passes:
            for oh, onode in other_nodes:
                dh = abs(self.node_height(onode)[0] - self.node_height(pss)[0])
                if node in self.crit_pts.peaks:
                    if pss not in snet[onode]:
                        snet.add_edge(onode, pss, dh=dh)
                        to_update.add(onode)
                        break
                elif node in self.crit_pts.pits:
                    if onode not in snet[pss]:
                        snet.add_edge(pss, onode, dh=dh)
                        to_update.add(onode)
                        break
        snet.remove_node(node)
        if node in self.crit_pts.pits:
            self.crit_pts.pits.remove(node)
        elif node in self.crit_pts.peaks:
            self.crit_pts.peaks.remove(node)
        if len(other_nodes) <= 1:
            snet.remove_node(c_pass)
            self.crit_pts.passes.remove(c_pass)
        for onode in to_update:
            assert node != onode
            assert c_pass != onode
            self.weight_peak_or_pit(onode, measure=measure)
        return to_update

    def get_peak_pit_region_area(self, node):
        region = self.get_minmax_region(node, self.crit_pts.passes)
        return len(region)

    def get_min_dh(self, node):
        snet = self.surf_network
        assert node in snet
        if node in self.crit_pts.pits: dd = snet.pred
        else: dd = snet.succ
        return min([dta['dh'] for dta in dd[node].values()])

    def get_max_dh(self, node):
        snet = self.surf_network
        assert node in snet
        if node in self.crit_pts.pits: dd = snet.pred
        else: dd = snet.succ
        return max([dta['dh'] for dta in dd[node].values()])

    def weight_peak_or_pit(self, node, measure=None):
        measure = measure or self.get_max_dh
        snet = self.surf_network
        snet.node[node] = measure(node)

    def weight_peaks_n_pits(self, measure=None):
        for pp in self.crit_pts.pits.union(self.crit_pts.peaks):
            self.weight_peak_or_pit(pp, measure=measure)
    
    def _get_surface_network(self):
        if self._surf_network:
            return self._surf_network
        passes = self.crit_pts.passes
        snet = netx.DiGraph()
        for p in passes:
            connected_peaks, connected_pits = self.nearest_extrema(p)
            pass_h = self.node_height(p)[0]
            for peak in connected_peaks:
                dh = self.node_height(peak)[0] - pass_h
                assert dh > 0.0
                snet.add_edge(peak, p, dh=dh)
            for pit in connected_pits:
                dh = pass_h - self.node_height(pit)[0]
                assert dh > 0.0
                snet.add_edge(p, pit, dh=dh)
        self._surf_network = snet
        self.regularize_surf_network(snet)
        return snet

    surf_network = property(_get_surface_network)

    def node_height(self, node, int sign=1):
        return (sign * self.arr[node], (sign * node[0], sign * node[1]))

    def get_reeb_graph(self):
        return get_reeb_graph(self.surf_network, self.crit_pts, self.node_height)

def is_pit_like(node, pits, passes, reeb, snet, node2h):
    if node in pits:
        return True
    elif node not in passes:
        # it's a peak
        return False
    reeb_nbrs = reeb._g[node]
    if len(reeb_nbrs) != 2:
        return False
    node_h = node2h(node)
    snet_nbrs = snet._g[node]
    for nbr in snet_nbrs:
        nbr_h = node2h(nbr)
        if nbr_h < node_h:
            return False
    return True

def is_peak_like(node, peaks, passes, reeb, snet, node2h):
    if node in peaks:
        return True
    elif node not in passes:
        # it's a pit
        return False
    reeb_nbrs = reeb._g[node]
    if len(reeb_nbrs) != 2:
        return False
    node_h = node2h(node)
    snet_nbrs = snet._g[node]
    for nbr in snet_nbrs:
        nbr_h = node2h(nbr)
        if nbr_h > node_h:
            return False
    return True

def is_unfinished(node, passes, gr):
    nbrs = gr._g[node]
    return node in passes and len(nbrs) != 3

def _is_peak_like(node, surf):
    return not surf.in_degree(node)

def _is_pit_like(node, surf):
    return not surf.out_degree(node)

def get_reeb_graph(surf_net, crit_pts, node2h):
    # make a copy of the surface network that we can destroy while creating
    # the reeb graph.
    surf = netx.DiGraph(surf_net)
    reeb = netx.Graph()
    peaks  = crit_pts.peaks
    passes = crit_pts.passes
    pits   = crit_pts.pits
    unfinished = []
    unfinished.extend(peaks)
    unfinished.extend(pits)
    unfinished.extend(passes)
    while unfinished:
        new_unfinished = []
        for node in unfinished:
            # if is_peak_like(node, peaks, passes, reeb, surf, node2h):
            if _is_peak_like(node, surf) and _is_pit_like(node, surf):
                continue
            if not (_is_peak_like(node, surf) or _is_pit_like(node, surf)):
                new_unfinished.append(node)
                continue
            elif _is_peak_like(node, surf):
                connected_passes = surf.succ[node]
                nearest_pass = max([(node2h(n), n) for n, dta in connected_passes.items()])[1]
            elif _is_pit_like(node, surf):
                connected_passes = surf.pred[node]
                nearest_pass = min([(node2h(n), n) for n, dta in connected_passes.items()])[1]
            other_passes = [n for n in connected_passes if n != nearest_pass]
            reeb.add_edge(node, nearest_pass)
            for other_pass in other_passes:
                # assert lower_pass in passes
                try:
                    surf.remove_edge(node, other_pass)
                except netx.NetworkXError:
                    surf.remove_edge(other_pass, node)
                dh = abs(node2h(nearest_pass)[0] - node2h(other_pass)[0])
                nearest_h = node2h(nearest_pass)
                other_h = node2h(other_pass)
                if nearest_h > other_h:
                    surf.add_edge(nearest_pass, other_pass, dh=dh)
                else:
                    surf.add_edge(other_pass, nearest_pass, dh=dh)
            try:
                surf.remove_edge(node, nearest_pass)
            except netx.NetworkXError:
                surf.remove_edge(nearest_pass, node)
        unfinished = new_unfinished
    return reeb

def _get_reeb_graph(surf_net, crit_pts, node2h):
    # make a copy of the surface network that we can destroy while creating
    # the reeb graph.
    ####
    # print "*** surface network ***", pformat(dict(surf_net._g))
    # print "crit pts: ", pformat(crit_pts)
    # for tp in crit_pts:
        # for node in crit_pts[tp]:
            # print "%s: %s," % (node, node2h(node))
    surf = surf_net.deepcopy()
    surf.order_by_key(keyfunc=node2h)
    reeb = graph()
    peaks  = crit_pts.peaks
    passes = crit_pts.passes
    pits   = crit_pts.pits
    unfinished = []
    unfinished.extend(peaks)
    unfinished.extend(passes)
    unfinished.extend(pits)
    while unfinished:
        new_unfinished = []
        for node in unfinished:
            if is_peak_like(node, peaks, passes, reeb, surf, node2h):
                # connect the peak to its highest pass.
                highest_pass = surf._g[node][-1]
                # assert highest_pass in passes
                reeb.add_edge(node, highest_pass)
                for lower_pass in surf._g[node][:-1]:
                    # assert lower_pass in passes
                    surf.remove_edge(node, lower_pass)
                    surf.add_edge(highest_pass, lower_pass)
                surf.remove_edge(node, highest_pass)
            elif is_pit_like(node, pits, passes, reeb, surf, node2h):
                # connect the pit to its lowest pass.
                lowest_pass = surf._g[node][0]
                # assert lowest_pass in passes
                reeb.add_edge(node, lowest_pass)
                for higher_pass in surf._g[node][1:]:
                    # assert higher_pass in passes
                    surf.remove_edge(node, higher_pass)
                    surf.add_edge(lowest_pass, higher_pass)
                surf.remove_edge(node, lowest_pass)
            elif is_unfinished(node, passes, reeb):
                new_unfinished.append(node)
        unfinished = new_unfinished
        # print 'unfinished: %s' % pformat(unfinished)
        # print 'reeb: %s' % pformat(dict(reeb._g))
        # print 'surf network: %s' % pformat(dict(surf._g))
        # res = raw_input('enter to continue, "q" to quit')
        # if res.lower() == 'q':
            # return reeb
    return reeb
