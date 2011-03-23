from collections import defaultdict

# from graph import Graph, DiGraph
from networkx import Graph, DiGraph

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

def join_split_tree_sparse(mesh, height_func, split=False):
    join_tree = DiGraph()
    # map of nodes to union-find set that it's in
    uf_map = {}
    # map from supernodes in the sparse j/s tree to all regular nodes in
    # supernode's arc.
    supernode_arc = {}
    # map of union-find set id to node to connect to in set.  When a join node
    # is created in the join_tree, connection_node_map keeps track of which
    # node should be connected in the join_tree.
    supernode_map = {}
    # a list of (height, node) tuples
    # XXX: store height_func() results in mesh?
    height_and_nodes = [(height_func(node), node) for node in mesh.nodes()]
    height_and_nodes.sort(reverse=not split)
    for h, n in height_and_nodes:
        uf_map[n] = set([n])
        supernode_map[id(uf_map[n])] = n
        connect_nbrs = []
        edges = []
        for nbr in mesh.neighbors_iter(n):
            # XXX: remove call to height_func()
            nbr_h = height_func(nbr)
            if nbr_h == h:
                raise ValueError("nbr_h == h!!!")
            if (split and nbr_h > h) or (not split and nbr_h < h) or (uf_map[nbr] is uf_map[n]):
                continue
            connection_nbr_uf = supernode_map[id(uf_map[nbr])]
            edges.append((connection_nbr_uf, n))
            connect_nbrs.append(nbr)
            uf_merge(uf_map, n, nbr)
        if len(edges) == 1:
            snode = supernode_map[id(uf_map[n])]
            supernode_arc[n] = snode
        elif len(edges) > 1:
            join_tree.add_edges_from(edges)
            for nbr in connect_nbrs:
                supernode_map[id(uf_map[nbr])] = n
    # add all the supernodes in the join_tree to the supernode_arcs.
    # do this before adding the global minimum.
    for node in join_tree.nodes():
        supernode_arc[node] = node
    # connect the global minimum to the join_tree's bottom
    tree_bottom = [n for (n,d) in join_tree.out_degree().items() if d == 0][0]
    global_min = height_and_nodes[-1][1]
    join_tree.add_edge(tree_bottom, global_min)
    return join_tree, supernode_arc

def splice_in_node(gr, start, newnode):
    succs = gr.successors(start)
    gr.add_edge(start, newnode)
    if succs:
        if len(succs) != 1:
            import pdb; pdb.set_trace()
        succ = succs[0]
        gr.add_edge(newnode, succ)
        gr.remove_edge(start, succ)

def get_regions_full(contour_tree):
    cpts = critical_points(contour_tree)
    peaks = cpts['peaks']
    pits = cpts['pits']
    passes = cpts['passes']
    cpts_flat = set(peaks.union(pits).union(passes))
    regions = {}
    for cpt in cpts_flat:
        region = [cpt]
        for lower_nbr in contour_tree.successors_iter(cpt):
            cur = lower_nbr
            while cur not in cpts_flat:
                region.append(cur)
                cur = contour_tree.successors(cur)[0]
            # add the last point
            region.append(cur)
            regions[(cpt, lower_nbr)] = region
    return regions

def get_regions_sparse(contour_tree, jnode2super, snode2super, height_func):
    jsuper2nodes = defaultdict(list)
    ssuper2nodes = defaultdict(list)
    for node in sorted(jnode2super, key=height_func):
        supern = jnode2super[node]
        jsuper2nodes[supern].append(node)
    for node in sorted(snode2super, key=height_func):
        supern = snode2super[node]
        ssuper2nodes[supern].append(node)
    regions = {}
    for higher_node in contour_tree:
        lower_nodes = contour_tree.successors(higher_node)
        for lower_node in lower_nodes:
            if higher_node in jsuper2nodes and lower_node in jsuper2nodes:
                # all info is in the jsuper2nodes map.
                region = jsuper2nodes[higher_node]
            elif higher_node in ssuper2nodes and lower_node in ssuper2nodes:
                # all info in ssuper2nodes map.
                region = ssuper2nodes[lower_node]
            elif higher_node in jsuper2nodes and lower_node in ssuper2nodes:
                # region is just the tail part of jsuper2nodes[higher_node]
                region = jsuper2nodes[higher_node]
                lower_h = height_func(lower_node)
                try:
                    idx = region.index(lower_node)
                except ValueError:
                    assert lower_h < height_func(region[0])
                    idx = 0
                region = region[idx:]
            elif higher_node in ssuper2nodes and lower_node in jsuper2nodes:
                # the complicated case...
                # find higher_node's superarc in jsuper2nodes
                hsuper = jnode2super[higher_node]
                hsuperarc = jsuper2nodes[hsuper]
                # and lower_node's superarc in ssuper2nodes
                lsuper = snode2super[lower_node]
                lsuperarc = ssuper2nodes[lsuper]
                # the region is the intersection between hsuprearc and lsuperarc
                region = sorted(set(hsuperarc).intersection(lsuperarc), key=height_func, reverse=True)
            regions[higher_node, lower_node] = region
    return regions

def rectify_join_split_trees(join, jarcs, split, sarcs, height_func):
    join_nodes = set(join.nodes())
    split_nodes = set(split.nodes())
    if join_nodes == split_nodes:
        return
    join_to_add = sorted(join_nodes.difference(split_nodes), key=height_func, reverse=True)
    split_to_add = sorted(split_nodes.difference(join_nodes), key=height_func)
    for jn in join_to_add:
        splice_in_node(split, sarcs[jn], jn)
    for sn in split_to_add:
        splice_in_node(join, jarcs[sn], sn)

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

def contour_tree_from_join_split(join, split, height_func):
    c_tree = DiGraph()
    leaves = join_split_peak_pit_nodes(join) + join_split_peak_pit_nodes(split)
    while len(leaves) > 1:
        leaf = leaves.pop()
        # if leaf == (301, 354):
            # import pdb; pdb.set_trace()
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

def sparse_contour_tree(mesh, height_func):
    join, jarcs = join_split_tree_sparse(mesh, height_func, split=False)
    split, sarcs = join_split_tree_sparse(mesh, height_func, split=True)
    rectify_join_split_trees(join, jarcs, split, sarcs, height_func)
    c_tree = contour_tree_from_join_split(join, split, height_func)
    return c_tree, get_regions_sparse(c_tree, jarcs, sarcs, height_func)

def contour_tree(mesh, height_func):
    join = join_split_tree(mesh, height_func, split=False)
    split = join_split_tree(mesh, height_func, split=True)
    return contour_tree_from_join_split(join, split, height_func)

def critical_points(ctree):
    crit_pts = {}
    in_deg = ctree.in_degree()
    out_deg = ctree.out_degree()
    # crit_pts['peaks']  = set([n for n in in_deg if in_deg[n] == 0 and out_deg[n] >= 1])
    # crit_pts['pits']   = set([n for n in out_deg if out_deg[n] == 0 and in_deg[n] >= 1])
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
    # return flattest(a, b, c, d)
    return envelope(a, b, c, d)
    # return AC

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
"""
def _join_split_tree(mesh, height_func, join_split_fac=1.0):
    join_tree = DiGraph()
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
            if nbr_h <= h or uf_map[nbr] is uf_map[n]:
                continue
            # uf_merge() can change the uf_map[nbr], and hence, the lowest
            # point in nbr's union-find set. So we add the edge first before doing the
            # union-find merge.
            lowest_in_nbr_uf = lowest_node_map[id(uf_map[nbr])]
            join_tree.add_edge(lowest_in_nbr_uf, n)
            uf_merge(uf_map, n, nbr)
            lowest_node_map[id(uf_map[nbr])] = n
    return join_tree
"""
