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

class graph:

    def __init__(self):
        from collections import defaultdict
        self._g = defaultdict(set)

    def add_edge(self, a, b):
        if b not in self._g[a]:
            self._g[a].add(b)
            self._g[b].add(a)

def mesh(arr):
    G = graph()
    nx, ny = arr.shape
    for i in range(nx):
        for j in range(ny):
            # connect along the cardinal directions
            G.add_edge((i,j), ((i+1)%nx, j))
            G.add_edge((i,j), ((i-1)%nx, j))
            G.add_edge((i,j), (i, (j+1)%ny))
            G.add_edge((i,j), (i, (j-1)%ny))
            a = arr[i,j]
            b = arr[i,(j+1)%ny]
            c = arr[(i+1)%nx, (j+1)%ny]
            d = arr[(i+1)%nx, j]
            if connect_diagonal(a, b, c, d) == 'ac':
                G.add_edge((i,j), ((i+1)%nx, (j+1)%ny))
            elif connect_diagonal(a, b, c, d) == 'bd':
                G.add_edge((i,(j+1)%ny), ((i+1)%nx, j))
            else:
                raise RuntimeError("invalid return value from connect_diagonal")
    return G

def order_neighbors(node, neighbors, shape):
    nx, ny = shape
    i, j = node
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
    assert set(ord_neighbors) == neighbors
    return ord_neighbors

def classify_nodes(arr, gr):
    classes = {'peaks': set(),
               'pits' : set(),
               'passes' : set(),
               }
    for node in gr._g:
        ordered_neighbors = order_neighbors(node, gr._g[node], arr.shape)
        diffs = [arr[n] - arr[node] for n in ordered_neighbors]
        assert any(diffs)
        diff_neg = diff_pos = 0.0
        for d in diffs:
            if d < 0:
                diff_neg += -d
            else:
                diff_pos += d
        n_change = 0
        for idx in range(len(diffs)):
            if diffs[idx] == 0.0 and diffs[idx-1] * diffs[(idx+1)%len(diffs)] < 0:
                n_change += 1
            elif diffs[idx-1] * diffs[idx] < 0:
                n_change += 1
        if n_change == 0:
            if diff_neg > 0 and diff_pos == 0.0:
                classes['peaks'].add(node)
            elif diff_pos > 0 and diff_neg == 0.0:
                classes['pits'].add(node)
        elif n_change == 4 and (diff_pos + diff_neg) > 0:
            classes['passes'].add(node)
    return classes

def sort_by_h(gr, arr):
    sgr = {}
    for node in gr._g:
        h_sort = sorted(gr._g[node], key=lambda n: arr[n])
        sgr[node] = h_sort
    return sgr

def surface_network(arr, gr, passes, peaks, pits):
    ctr_max = 4 * arr.shape[0]
    peaks_n_pits = peaks.union(pits)
    snet = graph()
    h_sorted_gr = sort_by_h(gr, arr)
    for p in passes:
        # get the two highest (lowest) neighbors to p.
        higher = h_sorted_gr[p][-2:]
        lower = h_sorted_gr[p][:2]
        nbrs_and_dir = ((higher, 'up'),
                        (lower , 'down'))
        for ns, dir in nbrs_and_dir:
            for n in ns:
                cur = n
                ctr = 0
                while cur not in peaks_n_pits and ctr < ctr_max:
                    if dir == 'up':
                        cur = h_sorted_gr[cur][-1]
                    else:
                        cur = h_sorted_gr[cur][0]
                    ctr += 1
                if ctr < ctr_max:
                    snet.add_edge(cur, p)
    return snet

def peak_pit_regions(arr, snet, classes):
    raise RuntimeError("not finished!!!")
    peaks = classes['peaks']
    pits = classes['pits']
    regions = {}
    for node in snet._g:
        if node in peaks:
            regions[node] = find_peak_region(node, snet, arr)
        elif node in pits:
            regions[node] = find_pit_region(node, snet, arr)
    return regions

def find_peak_region(node, snet, arr):
    raise RuntimeError("not finished!!!")
    for passes in snet._g[node]:
        h_sort = sorted(snet._g[node], key=lambda n: arr[n])
        highest_pass = h_sort[-1]
        h_val = arr[highest_pass]
        # while 
