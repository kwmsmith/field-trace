import networkx as nx

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

def join_split_tree(mesh, height_func):
    join_tree = nx.DiGraph()
    # map of nodes to union-find set that it's in
    uf_map = {}
    # map of union-find set id to lowest node in the set.
    lowest_node_map = {}
    # a list of (height, node) tuples
    height_and_nodes = [(height_func(node), node) for node in mesh.nodes()]
    height_and_nodes.sort(reverse=True)
    for h, n in height_and_nodes:
        uf_map[n] = set([n])
        lowest_node_map[id(uf_map[n])] = n
        for nbr in mesh.neighbors_iter(n):
            nbr_h = height_func(nbr)
            if nbr_h < h or uf_map[nbr] is uf_map[n]:
                continue
            # uf_merge() can change the uf_map[nbr] id value. Store it
            # for later retrieval.
            uf_map_nbr_id = id(uf_map[nbr])
            uf_merge(uf_map, n, nbr)
            join_tree.add_edge(lowest_node_map[uf_map_nbr_id], n)
            lowest_node_map[id(uf_map[nbr])] = n
    return join_tree

def peak_pit_nodes(join):
    in_deg = join.in_degree()
    return [n for n in in_deg if in_deg[n] == 0]

def pass_nodes(join):
    in_deg = join.in_degree()
    out_deg = join.out_degree()
    return [n for n in in_deg if in_deg[n] == 2 and out_deg[n] == 1]
