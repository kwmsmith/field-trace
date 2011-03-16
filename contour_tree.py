import networkx as _nx

def uf_merge(uf_map, key1, key2):
    s1 = uf_map[key1]
    s2 = uf_map[key2]
    if len(s1) <= len(s2):
        s2.update(s1)
        updated_set = s2
        keys_to_change = s1
    else:
        s1.update(s2)
        updated_set = s1
        keys_to_change = s2
    uf_map[key1] = uf_map[key2] = updated_set
    for key in keys_to_change:
        uf_map[key] = updated_set

def join_split_tree(mesh, height_func, join_split_fac=1.0):
    join_tree = _nx.DiGraph()
    # map of nodes to union-find set that it's in
    uf_map = {}
    # map of union-find set id to lowest node in the set.
    lowest_node_map = {}
    # a list of (height, node) tuples
    height_and_nodes = [(join_split_fac * height_func(node), node) for node in mesh.nodes()]
    height_and_nodes.sort(reverse=True)
    for h, n in height_and_nodes:
        uf_map[n] = set([n])
        lowest_node_map[id(uf_map[n])] = n
        for nbr in mesh.neighbors_iter(n):
            nbr_h = join_split_fac * height_func(nbr)
            if nbr_h < h or uf_map[nbr] is uf_map[n]:
                continue
            # uf_merge() can change the uf_map[nbr], and hence, the lowest
            # point in nbr's union-find set. So we add the edge first before doing the
            # union-find merge.
            lowest_in_nbr_uf = lowest_node_map[id(uf_map[nbr])]
            join_tree.add_edge(lowest_in_nbr_uf, n)
            uf_merge(uf_map, n, nbr)
            lowest_node_map[id(uf_map[nbr])] = n
    return join_tree

def join_split_peak_pit_nodes(join):
    in_deg = join.in_degree()
    return [n for n in in_deg if in_deg[n] == 0]

def join_split_pass_nodes(join):
    in_deg = join.in_degree()
    out_deg = join.out_degree()
    return [n for n in in_deg if in_deg[n] >= 2 and out_deg[n] == 1]

def is_leaf(node, join, split):
    return is_upper_leaf(node, join, split) or is_lower_leaf(node, join, split)

def is_upper_leaf(leaf, join, split):
    return not join.in_degree(leaf) and split.in_degree(leaf) == 1

def is_lower_leaf(leaf, join, split):
    return not split.in_degree(leaf) and join.in_degree(leaf) == 1

def reduce_graph(graph, node):
    preds = graph.predecessors(node)
    succs = graph.successors(node)
    if preds:
        assert len(preds) == 1
        if succs:
            assert len(succs) == 1
            graph.add_edge(preds[0], succs[0])
        graph.remove_edge(preds[0], node)
    elif succs:
        assert len(succs) == 1
        graph.remove_edge(node, succs[0])
    graph.remove_node(node)

def contour_tree(mesh, height_func):
    join = join_split_tree(mesh, height_func, join_split_fac=1.0)
    join.name = 'join'
    split = join_split_tree(mesh, height_func, join_split_fac=-1.0)
    split.name = 'split'
    c_tree = _nx.DiGraph()
    leaves = join_split_peak_pit_nodes(join) + join_split_peak_pit_nodes(split)
    while len(leaves) > 1:
        leaf = leaves.pop()
        if is_upper_leaf(leaf, join, split):
            this = join
            other = split
        else:
            this = split
            other = join
        this_succ = this.successors(leaf)
        assert len(this_succ) == 1
        nbr = this_succ[0]
        if height_func(leaf) > height_func(nbr):
            c_tree.add_edge(leaf, nbr)
        else:
            c_tree.add_edge(nbr, leaf)
        reduce_graph(this, leaf)
        reduce_graph(other, leaf)
        if is_leaf(nbr, join, split):
            leaves.append(nbr)
    return c_tree

def critical_points(ctree):
    crit_pts = {}
    in_deg = ctree.in_degree()
    out_deg = ctree.out_degree()
    crit_pts['peaks']  = set([n for n in in_deg if in_deg[n] == 0])
    crit_pts['pits']   = set([n for n in out_deg if out_deg[n] == 0])
    crit_pts['passes'] = set([n for n in out_deg if out_deg[n] >= 1 and in_deg[n] >= 2] +
                             [n for n in out_deg if out_deg[n] >= 2 and in_deg[n] >= 1])
    return crit_pts

AC = 0
BD = 1

def envelope(a, b, c, d):
    if b + d >= a + c:
        return BD
    else:
        return AC

def flattest(a, b, c, d):
    if (d - b)*(d - b) >= (a - c)*(a - c):
        return AC
    else:
        return BD

def connect_diagonal(a, b, c, d):
    """
    returns 'ac' or 'bd' indicating which nodes to connect in the square cell.
    
    """
    return flattest(a, b, c, d)
    # return envelope(a, b, c, d)
    # return AC

def make_mesh(arr):
    G = _nx.Graph()
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
            if connect_diagonal(a, b, c, d) == AC:
                G.add_edge((i,j), ((i+1)%nx, (j+1)%ny))
            elif connect_diagonal(a, b, c, d) == BD:
                G.add_edge((i,(j+1)%ny), ((i+1)%nx, j))
            else:
                raise RuntimeError("invalid return value from connect_diagonal")
    return G
