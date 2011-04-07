# cython: profile=True

from pprint import pprint, pformat

import networkx as netx

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

def _get_lowest_neighbor(h_sorted_mesh, node):
    return h_sorted_mesh[node][0]

def _get_highest_neighbor(h_sorted_mesh, node):
    return h_sorted_mesh[node][-1]

class TopoSurface(object):

    def __init__(self, arr):
        self.arr = arr
        self.mesh = self.get_mesh()
        self.h_sorted_mesh = sort_by_h(self.mesh, self.arr)
        self.crit_pts = self.get_crit_pts()
        self._surf_network = None

    def get_min_region(self, peak, passes):
        return self._get_minmax_region(peak, passes, sign=1)

    def get_max_region(self, peak, passes):
        return self._get_minmax_region(peak, passes, sign=-1)

    def _get_minmax_region(self, pit, passes, int sign):
        raise RuntimeError("test me!!!")
        import heapq
        cdef set region = set([pit])
        cdef list frontier = [(self.node_height(n, sign=sign), n) for n in self.mesh.neighbors(pit)]
        heapq.heapify(frontier)
        while True:
            h, n = heapq.heappop(frontier)
            region.add(n)
            if n in passes:
                break
            nbrs = self.mesh.neighbors(n)
            for nbr in nbrs:
                if nbr not in region:
                    heapq.heappush((self.node_height(nbr, sign=sign), nbr))
        return region

    def get_mesh(self):
        G = graph()
        nx, ny = self.arr.shape
        for i in range(nx):
            for j in range(ny):
                # connect along the cardinal directions
                G.add_edge((i,j), ((i+1)%nx, j))
                G.add_edge((i,j), ((i-1)%nx, j))
                G.add_edge((i,j), (i, (j+1)%ny))
                G.add_edge((i,j), (i, (j-1)%ny))
                a = self.arr[i,j]
                b = self.arr[i,(j+1)%ny]
                c = self.arr[(i+1)%nx, (j+1)%ny]
                d = self.arr[(i+1)%nx, j]
                if connect_diagonal(a, b, c, d) == AC:
                    G.add_edge((i,j), ((i+1)%nx, (j+1)%ny))
                elif connect_diagonal(a, b, c, d) == BD:
                    G.add_edge((i,(j+1)%ny), ((i+1)%nx, j))
                else:
                    raise RuntimeError("invalid return value from connect_diagonal")
        G.order_neighbors(self.arr.shape)
        return G

    def get_crit_pts(self):
        cdef double diff_neg, diff_pos, diff1, diff2
        cdef int n_change, idx
        from collections import namedtuple
        cps = namedtuple('crit_pts', 'peaks pits passes')
        crit_pts = {'peaks': set(),
                   'pits' : set(),
                   'passes' : set(),
                   }
        for node in self.mesh._g:
            nbrs = self.mesh._g[node]
            diffs = [self.arr[n] - self.arr[node] for n in nbrs]
            # assert any(diffs)
            diff_neg = diff_pos = 0.0
            for d in diffs:
                dd = d
                if dd < 0:
                    diff_neg += -dd
                else:
                    diff_pos += dd
            n_change = 0
            for idx in range(len(diffs)):
                diff1 = diffs[idx-1]
                diff2 = diffs[idx]
                # if diffs[idx] == 0.0 and diffs[idx-1] * diffs[(idx+1)%len(diffs)] < 0:
                    # n_change += 1
                # elif diffs[idx-1] * diffs[idx] < 0:
                if diff1 * diff2 < 0:
                    n_change += 1
            if n_change == 0:
                if diff_neg > 0 and diff_pos == 0.0:
                    crit_pts['peaks'].add(node)
                elif diff_pos > 0 and diff_neg == 0.0:
                    crit_pts['pits'].add(node)
            elif (n_change >= 4 and not n_change % 2) and (diff_pos + diff_neg) > 0:
                crit_pts['passes'].add(node)
        return cps(crit_pts['peaks'], crit_pts['pits'], crit_pts['passes'])
    
    def nearest_extrema(self, pss):
        pass_nbrs = self.mesh._g[pss]
        peaks_n_pits = self.crit_pts.peaks.union(self.crit_pts.pits)
        connected_peaks = set()
        connected_pits = set()
        for nbr in pass_nbrs:
            # go up from every neighbor.
            cur = nbr
            while True:
                if cur in self.crit_pts.peaks:
                    connected_peaks.add(cur)
                    break
                cur = _get_highest_neighbor(self.h_sorted_mesh, cur)
                if cur == pss or cur in pass_nbrs:
                    # returned to a previously visited spot; ignore these.
                    break
            # go down from every neighbor.
            cur = nbr
            while True:
                if cur in self.crit_pts.pits:
                    connected_pits.add(cur)
                    break
                cur = _get_lowest_neighbor(self.h_sorted_mesh, cur)
                if cur == pss or cur in pass_nbrs:
                    # returned to a previously visited spot; ignore these.
                    break
        return (connected_peaks, connected_pits)
    
    def get_surface_network(self):
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
        return snet

    surf_network = property(get_surface_network)

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
    print node, "node in passes: %s" % (node in passes), "len(nbrs) = %d" % len(nbrs)
    return node in passes and len(nbrs) != 3

def get_reeb_graph(surf_net, crit_pts, node2h):
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
