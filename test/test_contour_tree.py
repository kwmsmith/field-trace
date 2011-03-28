import networkx as _nx

from random import randint

import contour_tree as ct

from nose.tools import eq_, set_trace

from test_critical_point_network import random_periodic_upsample, visualize

def _set_nbr_height_equal(arr):
    nx, ny = arr.shape
    xeq = randint(1, nx-2)
    yeq = randint(1, ny-2)
    v1 = arr[xeq, yeq]
    v2 = arr[xeq, yeq+1]
    vavg = 0.5 * (v1+v2)
    arr[xeq, yeq] = arr[xeq, yeq+1] = vavg

class test_prune_regions(object):

    def setUp(self):
        self.mesh = ct.Graph()
        self.mesh.add_edges_from(mesh_edges)
        self.height_func = lambda n: n
        self.region_func = lambda n: len(n)

    def test_prune_regions(self):
        contour_tree = ct.contour_tree(self.mesh, self.height_func)
        ctree_copy = contour_tree.copy()
        ct.prune_regions(contour_tree, region_func=self.region_func, threshold=3, height_func=self.height_func)
        if 0:
            import pylab as pl
            pl.ion()
            _nx.draw(ctree_copy)
            pl.figure()
            _nx.draw(contour_tree)
            raw_input("enter to continue")
        set_trace()
        eq_(ctree_copy.adj, contour_tree.adj)

class test_contour_tree(object):

    def setUp(self):
        self.mesh = ct.Graph()
        self.mesh.add_edges_from(mesh_edges)
        self.height_func = lambda n: n

    def test_join_split_trees(self):
        eq_(len(self.mesh), 18)
        join = ct.join_split_tree(self.mesh, self.height_func)
        split = ct.join_split_tree(self.mesh, self.height_func, split=True)
        eq_(join.order(), split.order())
        eq_(sorted(join.nodes()), sorted(split.nodes()))
        eq_(sorted(ct.join_split_peak_pit_nodes(join)), [7, 8, 9, 10])
        eq_(sorted(ct.join_split_peak_pit_nodes(split)), [1, 2])
        eq_(sorted(ct.join_split_pass_nodes(join)), [4, 5, 6])
        eq_(sorted(ct.join_split_pass_nodes(split)), [3])

        compare_adj(join.adj, join_adj)
        compare_adj(split.adj, split_adj)

        if 0:
            import pylab as pl
            pl.ion()
            nx_join = join.to_networkx()
            _nx.draw_shell(nx_join)
            pl.title('join')
            pl.figure()
            nx_split = split.to_networkx()
            _nx.draw_shell(nx_split)
            pl.title('split')
            raw_input("enter to continue")

    def test_contour_tree(self):
        join = ct.join_split_tree(self.mesh, self.height_func)
        split = ct.join_split_tree(self.mesh, self.height_func, split=True)
        contour_tree = ct.contour_tree(self.mesh, self.height_func)
        crit_pts = ct.critical_points(contour_tree)
        eq_(sorted(crit_pts['peaks']), sorted(ct.join_split_peak_pit_nodes(join)))
        eq_(sorted(crit_pts['pits']), sorted(ct.join_split_peak_pit_nodes(split)))
        eq_(sorted(crit_pts['passes']), sorted([3, 4, 5, 6]))
        regions, sparse_tree = ct.get_regions_full(contour_tree, sparse_tree=True)

        if 0:
            import pylab as pl
            pl.ion()
            try:
                _nx.draw(sparse_tree.to_networkx())
            except AttributeError:
                _nx.draw(sparse_tree)
            pl.figure()
            try:
                _nx.draw(contour_tree.to_networkx())
            except AttributeError:
                _nx.draw(contour_tree)
            raw_input("enter to continue")

    def _test_contour_tree_sparse(self):
        join, join_arcs = ct.join_split_tree_sparse(self.mesh, self.height_func)
        eq_(set(join_arcs.keys()), set(self.mesh.nodes()))
        split, split_arcs = ct.join_split_tree_sparse(self.mesh, self.height_func, split=True)
        eq_(set(split_arcs.keys()), set(self.mesh.nodes()))
        ct.rectify_join_split_trees(join, join_arcs, split, split_arcs, self.height_func)
        eq_(sorted(join.nodes()), sorted(split.nodes()))
        contour_tree, regions = ct.sparse_contour_tree(self.mesh, self.height_func)
        crit_pts = ct.critical_points(contour_tree)
        eq_(sorted(crit_pts['peaks']), sorted(ct.join_split_peak_pit_nodes(join)))
        eq_(sorted(crit_pts['pits']), sorted(ct.join_split_peak_pit_nodes(split)))
        eq_(sorted(crit_pts['passes']), sorted([3, 4, 5, 6]))
        if 0:
            import pylab as pl
            pl.ion()
            nx_join = join.to_nx()
            _nx.draw(nx_join)
            pl.title('join')
            pl.figure()
            nx_split = split.to_nx()
            _nx.draw(nx_split)
            pl.title('split')
            pl.figure()
            nx_ctree = contour_tree.to_nx()
            _nx.draw(nx_ctree)
            raw_input("enter to continue")
        eq_(regions,
                {
                    (6, 4): [4.9, 6],
                    (3, 2): [2],
                    (7, 6): [7],
                    (3, 1): [1, 2.1],
                    (5, 4): [4.6, 5],
                    (10, 5): [6.1, 7.2, 8.3, 10],
                    (4, 3): [3, 4],
                    (9, 5): [6.5, 9],
                    (8, 6): [6.9, 8],
                    })

    NN = 32
    def test_arr_full(self):
        arr = random_periodic_upsample(self.NN, 4, seed=1)
        for _ in range(4):
            _set_nbr_height_equal(arr)
        mesh = ct.make_mesh(arr)
        def height_func(n):
            return (arr[n], n)
        contour_tree = ct.contour_tree(mesh, height_func)
        cpts = ct.critical_points(contour_tree)
        peaks = cpts['peaks']
        passes = cpts['passes']
        pits = cpts['pits']
        print "peaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        print "len(crit_pts) = %d" % (len(peaks) + len(pits) + len(passes))
        print 'tot points covered: %d' % len(contour_tree)
        coverage = set()
        regions = [D['arc'] for (n1, n2, D) in contour_tree.edges_iter(data=True)]
        for r in regions:
            coverage.update(r)
        eq_(len(coverage), arr.size)
        regions = ct.get_regions_full(contour_tree)
        if 0:
            vis(arr, height_func=height_func, crit_pts=cpts, regions=regions)

def vis(arr, height_func, crit_pts, regions, step=True, new_fig=False):
    filled_arr = arr.copy()
    def hf((a,b)): return height_func(a), height_func(b)
    for region in sorted(regions, key=hf, reverse=True):
        # filled_arr = arr.copy()
        X = [_[0] for _ in regions[region]]
        Y = [_[1] for _ in regions[region]]
        filled_arr[X,Y] = 2*arr.max()
        if step:
            visualize(filled_arr, crit_pts=crit_pts, ncontours=None, cmap='gray', new_fig=new_fig)
            raw_input("enter to continue")
    if not step:
        visualize(filled_arr, crit_pts=crit_pts, ncontours=None, cmap='gray', new_fig=True)
        raw_input("enter to continue")


mesh_edges = [
        (1, 2.1),
        (1, 9),
        (1, 8),

        (2, 3),
        (2, 4.6),
        (2, 5),
        (2, 4),
        (2, 6),
        (2, 4.9),

        (2.1, 9),
        (2.1, 3),
        (2.1, 8),

        (3, 9),
        (3, 6.5),
        (3, 5),
        (3, 4.9),
        (3, 6.9),
        (3, 8),

        (4, 4.6),
        (4, 6.1),
        (4, 10),
        (4, 7),
        (4, 6),
        (4, 5),

        (4.6, 5),

        (4.9, 6.9),
        (4.9, 6),

        (5, 6.5),
        (5, 9),
        (5, 10),
        (5, 7.2),

        (6, 8),
        (6, 6.9),
        (6, 7),

        (6.1, 8.3),
        (6.1, 10),
        (6.1, 7.2),

        (6.5, 9),

        (6.9, 8),

        (7.2, 10),
        (7.2, 8.3),

        (8.3, 10),
    ]

def compare_adj(adj1, adj2):
    for x in adj1:
        eq_(set(adj1[x]), set(adj2[x]))

join_adj = {
        1: set([]),
        2: set([1]),
        2.1: set([2]),
        3: set([2.1]),
        4: set([3]),
        4.6: set([4]),
        4.9: set([4]),
        5: set([4.6]),
        6: set([4.9]),
        6.1: set([5]),
        6.5: set([5]),
        6.9: set([6]),
        7: set([6]),
        7.2: set([6.1]),
        8: set([6.9]),
        8.3: set([7.2]),
        9: set([6.5]),
        10: set([8.3]),
        }

split_adj = {
    1: set([2.1]),
    2: set([3]),
    2.1: set([3]),
    3: set([4]),
    4: set([4.6]),
    4.6: set([4.9]),
    4.9: set([5]),
    5: set([6]),
    6: set([6.1]),
    6.1: set([6.5]),
    6.5: set([6.9]),
    6.9: set([7]),
    7: set([7.2]),
    7.2: set([8]),
    8: set([8.3]),
    8.3: set([9]),
    9: set([10]),
    10: set([]),
    }
