import networkx as _nx

from random import randint

import contour_tree as ct

from nose.tools import eq_, ok_, set_trace

from test_critical_point_network import random_periodic_upsample, visualize

import pylab as pl
pl.ion()

def _set_nbr_height_equal(arr):
    nx, ny = arr.shape
    xeq = randint(1, nx-2)
    yeq = randint(1, ny-2)
    v1 = arr[xeq, yeq]
    v2 = arr[xeq, yeq+1]
    vavg = 0.5 * (v1+v2)
    arr[xeq, yeq] = arr[xeq, yeq+1] = vavg

class test_hierarchical_cluster(object):
    
    def setUp(self):
        self.seed_pts = [(0,0), (1,1), (0,1), (511, 511),
                         (255, 255), (255, 256),
                         (511, 0), (510, 3),
                         (0, 511),
                         (511, 256), (0, 256),
                         (256, 511), (256, 0),
                         ]
        self.nx = self.ny = 512

    def test_hc(self):
        Z = ct.hierarchical_cluster(self.seed_pts, self.nx, self.ny)
        clusters = ct.fclusters(Z, self.seed_pts, nclusters=4)
        eq_(clusters,
                [
                set([(255, 256), (255, 255)]),
                set([(256, 0), (256, 511)]),
                set([(0, 256), (511, 256)]),
                set([(0, 1), (0, 0), (510, 3), (511, 511), (0, 511), (511, 0), (1, 1)]),
                ])

    def test_bounding_box(self):
        Z = ct.hierarchical_cluster(self.seed_pts, self.nx, self.ny)
        clusters = ct.fclusters(Z, self.seed_pts, nclusters=4)
        bbs = []
        for cluster in clusters:
            bbs.append(ct.bounding_box(cluster, self.nx, self.ny))

class test_prune_regions(object):

    def setUp(self):
        self.mesh = ct.Graph()
        self.mesh.add_edges_from(mesh_edges)
        self.height_func = lambda n: n
        self.region_func = lambda n: len(n)

    def test_prune_regions(self):
        contour_tree = ct.contour_tree(self.mesh, self.height_func)
        crit_pts0 = ct.critical_points(contour_tree)
        pits, passes, peaks = crit_pts0.pits, crit_pts0.passes, crit_pts0.peaks
        regions = ct.get_regions(contour_tree)
        domain = set()
        for r in regions.values():
            domain.update(r)
        ctree_copy = contour_tree.copy()
        thresh = 0
        while len(contour_tree) > 2:
            ct.prune_regions(contour_tree, region_func=self.region_func, threshold=thresh)
            regions = ct.get_regions(contour_tree)
            domain_pruned = set()
            for r in regions.values():
                domain_pruned.update(r)
            pruned_crit_pts = ct.critical_points(contour_tree)
            ppeaks = pruned_crit_pts.peaks
            ppits = pruned_crit_pts.pits
            ppasses = pruned_crit_pts.passes
            eq_(domain, domain_pruned)
            eq_(len(ppeaks)+len(ppits)-len(ppasses), len(pits)+len(peaks)-len(passes))
            ok_(ppeaks <= peaks)
            ok_(ppasses <= passes)
            ok_(ppits <= pits)
            if 0:
                import pylab as pl
                pl.ion()
                _nx.draw(ctree_copy)
                pl.figure()
                _nx.draw(contour_tree)
                raw_input("enter to continue")
            thresh += 1

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
        eq_(sorted(crit_pts.peaks), sorted(ct.join_split_peak_pit_nodes(join)))
        eq_(sorted(crit_pts.pits), sorted(ct.join_split_peak_pit_nodes(split)))
        eq_(sorted(crit_pts.passes), sorted([3, 4, 5, 6]))

        if 0:
            import pylab as pl
            pl.ion()
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
        eq_(sorted(crit_pts.peaks), sorted(ct.join_split_peak_pit_nodes(join)))
        eq_(sorted(crit_pts.pits), sorted(ct.join_split_peak_pit_nodes(split)))
        eq_(sorted(crit_pts.passes), sorted([3, 4, 5, 6]))
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
        def region_func(r):
            hs = [height_func(p) for p in r]
            return max(hs)[0] - min(hs)[0]
        def region_area_func(r):
            return len(r)
        contour_tree = ct.contour_tree(mesh, height_func)
        def remove_edge_cb(region, interior, leaf):
            h = height_func(interior)[0]
            for pt in region:
                arr[pt] = h
        ct.prune_regions(contour_tree,
                region_func=region_area_func,
                # threshold=(arr.max()-arr.min())/4.0,
                threshold=3,
                remove_edge_cb=remove_edge_cb)
        cpts = ct.critical_points(contour_tree)
        peaks = cpts.peaks
        passes = cpts.passes
        pits = cpts.pits
        print "peaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        print "len(crit_pts) = %d" % (len(peaks) + len(pits) + len(passes))
        print 'tot points covered: %d' % len(contour_tree)
        coverage = set()
        regions = ct.get_regions(contour_tree)
        for r in regions.values():
            coverage.update(r)
        eq_(len(coverage), arr.size)
        if 0:
            vis(arr, height_func=height_func, crit_pts=cpts, regions=regions)

def vis(arr, height_func, crit_pts, regions, step=True, new_fig=False):
    filled_arr = arr.copy()
    def hf((a,b)): return height_func(a), height_func(b)
    for region in sorted(regions, key=hf, reverse=True):
        if step:
            visualize(filled_arr, crit_pts=crit_pts, ncontours=None, cmap='gray', new_fig=new_fig)
            raw_input("enter to continue")
        X = [_[0] for _ in regions[region]]
        Y = [_[1] for _ in regions[region]]
        filled_arr[X,Y] = 2*arr.max()
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
