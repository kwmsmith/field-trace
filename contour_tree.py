from _graph import Graph, DiGraph
# from networkx import Graph, DiGraph

def contour_tree(mesh, height_func):
    join = join_split_tree(mesh, height_func, split=False)
    split = join_split_tree(mesh, height_func, split=True)
    c_tree = contour_tree_from_join_split(join, split, height_func)
    supertree = augmented_tree_to_supertree(c_tree)
    return supertree

def join_split_tree(mesh, height_func, split=False):
    join_tree = DiGraph()
    # map of nodes to union-find set that it's in
    uf_map = {}
    # map of union-find set id to lowest node in the set.
    lowest_node_map = {}
    # a list of (height, node) tuples
    height_and_nodes = [(height_func(node), node) for node in mesh.nodes()]
    height_and_nodes.sort(reverse=not split)
    for h, n in height_and_nodes:
        uf_map[n] = set([n])
        lowest_node_map[id(uf_map[n])] = n
        for nbr in mesh.neighbors_iter(n):
            nbr_h = height_func(nbr)
            if nbr_h == h:
                raise ValueError("nbr_h == h!!!")
            if (split and nbr_h > h) or (not split and nbr_h < h) or (uf_map[nbr] is uf_map[n]):
                continue
            # uf_merge() can change the uf_map[nbr], and hence, the lowest
            # point in nbr's union-find set. So we add the edge first before doing the
            # union-find merge.
            lowest_in_nbr_uf = lowest_node_map[id(uf_map[nbr])]
            join_tree.add_edge(lowest_in_nbr_uf, n)
            uf_merge(uf_map, n, nbr)
            lowest_node_map[id(uf_map[nbr])] = n
    return join_tree

def contour_tree_from_join_split(join, split, height_func):
    c_tree = DiGraph()
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

def augmented_tree_to_supertree(contour_tree):
    '''
    takes a fully augmented contour tree and returns an unaugmented contour
    tree with the edge arcs as edge attributes.

    Note: always returns a networkx Graph.
    '''
    import networkx as _nx
    cpts = critical_points(contour_tree)
    peaks = cpts.peaks
    pits = cpts.pits
    passes = cpts.passes
    cpts_flat = set(peaks.union(pits).union(passes))
    edge2regions = []
    for cpt in cpts_flat:
        for lower_nbr in contour_tree.successors(cpt):
            region = set([cpt])
            cur = lower_nbr
            # if not isinstance(cpt, tuple):
                # import pdb; pdb.set_trace()
            while cur not in cpts_flat:
                region.add(cur)
                cur = contour_tree.successors(cur)[0]
            # add the last point
            region.add(cur)
            edge2regions.append((cpt, cur, dict(arc=region)))
    supertree = _nx.DiGraph(edge2regions)
    return supertree

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

def splice_in_node(gr, start, newnode):
    succs = gr.successors(start)
    gr.add_edge(start, newnode)
    if succs:
        if len(succs) != 1:
            raise ValueError("len(succs) != 1: %s" % (succs,))
        succ = succs[0]
        gr.add_edge(newnode, succ)
        gr.remove_edge(start, succ)

def get_regions(contour_tree):
    regions = dict([((n1, n2), D['arc']) for (n1, n2, D) in contour_tree.edges_iter(data=True)])
    return regions

def join_split_peak_pit_nodes(join):
    in_deg = join.in_degree()
    return [n for n in in_deg if in_deg[n] == 0]

def join_split_pass_nodes(join):
    in_deg = join.in_degree()
    out_deg = join.out_degree()
    return [n for n in in_deg if in_deg[n] >= 2 and out_deg[n] == 1]

def is_leaf(node, join, split):
    return is_upper_leaf(node, join, split) or is_lower_leaf(node, join, split)

def is_upper_leaf(node, join, split):
    return not join.in_degree(node) and split.in_degree(node) == 1

def is_lower_leaf(node, join, split):
    return not split.in_degree(node) and join.in_degree(node) == 1

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

def critical_points(ctree):
    from collections import namedtuple
    crit_pts = namedtuple('crit_pts', 'peaks passes pits')
    in_deg = ctree.in_degree()
    out_deg = ctree.out_degree()
    return crit_pts(
            peaks = set([n for n in in_deg if in_deg[n] == 0]),
            pits = set([n for n in out_deg if out_deg[n] == 0]),
            passes = set([n for n in out_deg if (out_deg[n] + in_deg[n] > 2)]),
            )

def _critical_points(ctree):
    # this version is slower but could be made faster by removing
    # ctree.in_degree() and ctree.out_degree() overhead
    deg = ctree.degree()
    cpts = [n for n in deg if deg[n] != 2]
    crit_pts_dict = {}
    in_deg = dict([(cp, ctree.in_degree(cp)) for cp in cpts])
    out_deg = dict([(cp, ctree.out_degree(cp)) for cp in cpts])
    peaks = crit_pts_dict['peaks'] = set([n for n in in_deg if in_deg[n] == 0])
    pits = crit_pts_dict['pits'] = set([n for n in out_deg if out_deg[n] == 0])
    crit_pts_dict['passes'] = set(cpts).difference(peaks.union(pits))
    return crit_pts_dict

def make_mesh(arr):
    G = Graph()
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
            con_dp = connect_diagonal(a, b, c, d)
            if con_dp == AC:
                G.add_edge((i,j), ((i+1)%nx, (j+1)%ny))
            elif con_dp == BD:
                G.add_edge((i,(j+1)%ny), ((i+1)%nx, j))
            else:
                raise RuntimeError("invalid return value from connect_diagonal")
    return G

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
    # return flattest(a, b, c, d)
    return envelope(a, b, c, d)
    # return AC

def get_leaf_edges(c_tree, data=True):
    in_deg = c_tree.in_degree()
    out_deg = c_tree.out_degree()
    edges = []
    for e in c_tree.edges_iter(data=data):
        if data:
            n1, n2, d = e
        else:
            n1, n2 = e
        if in_deg[n1] == 0:
            edges.append(e)
        if out_deg[n2] == 0:
            edges.append(e)
    return edges

def interior_exterior(c_tree, edge):
    n1, n2 = edge
    if c_tree.in_degree(n1) == 0 or c_tree.out_degree(n1) == 0:
        # n1 is exterior
        # n2 is interior
        return (n2, n1)
    elif c_tree.in_degree(n2) == 0 or c_tree.out_degree(n2) == 0:
        # flip
        return (n1, n2)
    else:
        raise ValueError("can't figure it out!")

def is_upper_edge(c_tree, edge):
    return c_tree.in_degree(edge[0]) == 0

def is_lower_edge(c_tree, edge):
    return c_tree.out_degree(edge[1]) == 0

def is_last_up_down(c_tree, edge, interior_node):
    # can be sped up, I'm sure...
    if ((c_tree.in_degree(interior_node) == 1  and is_upper_edge(c_tree, edge)) or 
        (c_tree.out_degree(interior_node) == 1 and is_lower_edge(c_tree, edge))):
        return True
    return False

def pred_edges(c_tree, node):
    pes = set([(n, node) for n in c_tree.pred[node]])
    return pes

def succ_edges(c_tree, node):
    ses = set([(node, n) for n in c_tree.succ[node]])
    return ses

def get_edge_region(c_tree, edge):
    return c_tree.edge[edge[0]][edge[1]]['arc']

def set_edge_region(c_tree, edge, region):
    c_tree.edge[edge[0]][edge[1]]['arc'] = region

def remove_edge_merge_regions(c_tree, leaf_edge, interior, leaf, already_collapsed, cb=None):
    p_edges = pred_edges(c_tree, interior)
    s_edges = succ_edges(c_tree, interior)
    if leaf_edge in p_edges:
        p_edges.remove(leaf_edge)
        sibling = p_edges.pop()
    elif leaf_edge in s_edges:
        s_edges.remove(leaf_edge)
        sibling = s_edges.pop()
    leaf_region = get_edge_region(c_tree, leaf_edge)
    sib_region = get_edge_region(c_tree, sibling)
    c_tree.remove_edge(*leaf_edge)
    c_tree.remove_node(leaf)
    set_edge_region(c_tree, sibling, sib_region.union(leaf_region))
    already_collapsed.add(leaf_edge)

def is_regular_node(c_tree, node):
    return c_tree.in_degree(node) == 1 and c_tree.out_degree(node) == 1

def node_collapse_merge_regions(c_tree, interior, already_collapsed):
    p_edges = pred_edges(c_tree, interior)
    s_edges = succ_edges(c_tree, interior)
    assert len(p_edges) == len(s_edges) == 1
    p_edge = p_edges.pop()
    s_edge = s_edges.pop()
    assert p_edge[1] == s_edge[0] == interior
    p_region = get_edge_region(c_tree, p_edge)
    s_region = get_edge_region(c_tree, s_edge)
    new_region = p_region.union(s_region)
    new_edge = p_edge[0], s_edge[1]
    c_tree.remove_node(interior)
    c_tree.add_edge(*new_edge, arc=new_region)
    already_collapsed.add(p_edge)
    already_collapsed.add(s_edge)
    return new_edge

def prune_regions(c_tree, region_func, threshold, height_func):
    '''
    prunes all edges in the supertree c_tree that have
    region_func(edge region) <= threshold.

    '''
    import heapq
    # load up queue
    already_collapsed = set()
    collapse_record = []
    leaf_edges = get_leaf_edges(c_tree)
    leaf_edges = [(region_func(D['arc']), (n1, n2)) for (n1, n2, D) in leaf_edges]
    heapq.heapify(leaf_edges)
    while leaf_edges and len(c_tree) > 2:
        priority, leaf_edge = heapq.heappop(leaf_edges)
        if priority > threshold:
            return
        if leaf_edge in already_collapsed:
            continue
        interior, leaf = interior_exterior(c_tree, leaf_edge)
        if is_last_up_down(c_tree, leaf_edge, interior):
            continue
        remove_edge_merge_regions(c_tree, leaf_edge, interior, leaf, already_collapsed)
        collapse_record.append(leaf_edge)
        if is_regular_node(c_tree, interior):
            new_leaf_edge = node_collapse_merge_regions(c_tree, interior, already_collapsed)
            if is_upper_edge(c_tree, new_leaf_edge) or is_lower_edge(c_tree, new_leaf_edge):
                u,v = new_leaf_edge
                region = get_edge_region(c_tree, new_leaf_edge)
                heapq.heappush(leaf_edges, (region_func(region), new_leaf_edge))
