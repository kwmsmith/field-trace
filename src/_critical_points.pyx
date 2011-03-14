def envelope(a, b, c, d):
    if b + d >= a + c:
        return 'bd'
    else:
        return 'ac'

def flattest(a, b, c, d):
    if (d - b)**2 >= (a - c)**2:
        return 'ac'
    else:
        return 'bd'

def connect_diagonal(a, b, c, d):
    """
    returns 'ac' or 'bd' indicating which nodes to connect in the square cell.
    
    """
    return flattest(a, b, c, d)
    # return envelope(a, b, c, d)
    # return 'ac'

class graph(object):

    def __init__(self):
        from collections import defaultdict
        self._g = defaultdict(list)

    def add_edge(self, a, b):
        if b not in self._g[a]:
            self._g[a].append(b)
            self._g[b].append(a)

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

def sort_by_h(gr, arr):
    sgr = {}
    for node in gr._g:
        h_sort = sorted(gr._g[node], key=lambda n: arr[n])
        sgr[node] = h_sort
    return sgr


class TopoSurface(object):

    def __init__(self, arr):
        self.arr = arr
        self.mesh = self.get_mesh()
        self.crit_pts = self.get_crit_pts()
        self.surf_network = self.get_surface_network()

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
                if connect_diagonal(a, b, c, d) == 'ac':
                    G.add_edge((i,j), ((i+1)%nx, (j+1)%ny))
                elif connect_diagonal(a, b, c, d) == 'bd':
                    G.add_edge((i,(j+1)%ny), ((i+1)%nx, j))
                else:
                    raise RuntimeError("invalid return value from connect_diagonal")
        G.order_neighbors(self.arr.shape)
        return G

    def get_crit_pts(self):
        crit_pts = {'peaks': set(),
                   'pits' : set(),
                   'passes' : set(),
                   }
        for node in self.mesh._g:
            nbrs = self.mesh._g[node]
            diffs = [self.arr[n] - self.arr[node] for n in nbrs]
            assert any(diffs)
            diff_neg = diff_pos = 0.0
            for d in diffs:
                if d < 0:
                    diff_neg += -d
                else:
                    diff_pos += d
            n_change = 0
            for idx in range(len(diffs)):
                # if diffs[idx] == 0.0 and diffs[idx-1] * diffs[(idx+1)%len(diffs)] < 0:
                    # n_change += 1
                # elif diffs[idx-1] * diffs[idx] < 0:
                if diffs[idx-1] * diffs[idx] < 0:
                    n_change += 1
            if n_change == 0:
                if diff_neg > 0 and diff_pos == 0.0:
                    crit_pts['peaks'].add(node)
                elif diff_pos > 0 and diff_neg == 0.0:
                    crit_pts['pits'].add(node)
            elif n_change == 4 and (diff_pos + diff_neg) > 0:
                crit_pts['passes'].add(node)
        return crit_pts

    def get_surface_network(self):
        peaks = self.crit_pts['peaks']
        pits = self.crit_pts['pits']
        passes = self.crit_pts['passes']
        ctr_max = 4 * self.arr.shape[0]
        peaks_n_pits = peaks.union(pits)
        snet = graph()
        h_sorted_mesh = sort_by_h(self.mesh, self.arr)
        for p in passes:
            # get the two highest (lowest) neighbors to p.
            higher, lower = self.separated_pass_nbrs(p)
            nbrs_and_dir = ((higher, 'up'),
                            (lower , 'down'))
            for ns, dir in nbrs_and_dir:
                for n in ns:
                    cur = n
                    ctr = 0
                    while cur not in peaks_n_pits and ctr < ctr_max:
                        if dir == 'up':
                            cur = h_sorted_mesh[cur][-1]
                        else:
                            cur = h_sorted_mesh[cur][0]
                        ctr += 1
                    if ctr < ctr_max:
                        snet.add_edge(cur, p)
        return snet

    def separated_pass_nbrs(self, pss):
        """
        returns neighbors of `pss` that sample different min & max regions;
        regions are all separated.

        """
        nbrs = self.mesh._g[pss]
        higher = []
        lower = []
        diffs = [self.arr[n] - self.arr[pss] for n in nbrs]
        for idx in range(len(diffs)):
            if diffs[idx-1] > 0 and diffs[idx] < 0:
                higher.append(nbrs[idx-1])
            elif diffs[idx-1] < 0 and diffs[idx] > 0:
                lower.append(nbrs[idx-1])
            elif diffs[idx] == 0.0:
                raise RuntimeError("gotcha: %s" % diffs)
        return higher, lower
