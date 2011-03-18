from pprint import pprint, pformat
import networkx as nx

import contour_tree as ct

from nose.tools import eq_, ok_, set_trace

from test_critical_point_network import random_periodic_upsample, visualize

class test_contour_tree(object):

    def setUp(self):
        self.mesh = nx.Graph()
        self.mesh.add_edges_from(mesh_edges)
        self.height_func = lambda n: n

    def test_join_split_trees(self):
        eq_(len(self.mesh), 18)
        join = ct.join_split_tree(self.mesh, self.height_func)
        split = ct.join_split_tree(self.mesh, self.height_func, join_split_fac=-1.0)
        eq_(join.order(), split.order())
        eq_(sorted(join.nodes()), sorted(split.nodes()))
        eq_(sorted(ct.join_split_peak_pit_nodes(join)), [7, 8, 9, 10])
        eq_(sorted(ct.join_split_peak_pit_nodes(split)), [1, 2])
        eq_(sorted(ct.join_split_pass_nodes(join)), [4, 5, 6])
        eq_(sorted(ct.join_split_pass_nodes(split)), [3])

        eq_(join.adj, join_adj)

        eq_(split.adj, split_adj)

        if 0:
            import pylab as pl
            pl.ion()
            nx.draw_shell(join)
            pl.title('join')
            pl.figure()
            nx.draw_shell(split)
            pl.title('split')
            raw_input("enter to continue")

    def test_contour_tree(self):
        join = ct.join_split_tree(self.mesh, self.height_func)
        split = ct.join_split_tree(self.mesh, self.height_func, join_split_fac=-1.0)
        contour_tree = ct.contour_tree(self.mesh, self.height_func)
        crit_pts = ct.critical_points(contour_tree)
        eq_(sorted(crit_pts['peaks']), sorted(ct.join_split_peak_pit_nodes(join)))
        eq_(sorted(crit_pts['pits']), sorted(ct.join_split_peak_pit_nodes(split)))
        eq_(sorted(crit_pts['passes']), sorted([3, 4, 5, 6]))

    def test_contour_tree_sparse(self):
        join, join_arcs = ct.join_split_tree_sparse(self.mesh, self.height_func)
        eq_(set(join_arcs.keys()), set(self.mesh.nodes()))
        split, split_arcs = ct.join_split_tree_sparse(self.mesh, self.height_func, join_split_fac=-1.0)
        eq_(set(split_arcs.keys()), set(self.mesh.nodes()))
        ct.rectify_join_split_trees(join, join_arcs, split, split_arcs, self.height_func)
        eq_(sorted(join.nodes()), sorted(split.nodes()))
        contour_tree, regions = ct.contour_tree(self.mesh, self.height_func, sparse=True)
        crit_pts = ct.critical_points(contour_tree)
        eq_(sorted(crit_pts['peaks']), sorted(ct.join_split_peak_pit_nodes(join)))
        eq_(sorted(crit_pts['pits']), sorted(ct.join_split_peak_pit_nodes(split)))
        eq_(sorted(crit_pts['passes']), sorted([3, 4, 5, 6]))
        if 0:
            import pylab as pl
            pl.ion()
            nx.draw(join)
            pl.title('join')
            pl.figure()
            nx.draw(split)
            pl.title('split')
            pl.figure()
            nx.draw(contour_tree)
            raw_input("enter to continue")
        eq_(regions,
                {
                    1: [1, 2.1],
                    2: [2],
                    3: [3],
                    4: [4],
                    5: [5, 4.6],
                    6: [6, 4.9],
                    7: [7],
                    8: [8, 6.9],
                    9: [9, 6.5],
                    10: [10, 8.3, 7.2, 6.1],
                    })

    def test_arr(self):
        arr = random_periodic_upsample(64, 8, seed=None)
        mesh = ct.make_mesh(arr)
        def height_func(n):
            return arr[n]
        contour_tree, regions = ct.contour_tree(mesh, height_func, sparse=True)
        cpts = ct.critical_points(contour_tree)
        peaks = cpts['peaks']
        passes = cpts['passes']
        pits = cpts['pits']
        print arr[0,0]
        print "peaks + pits - passes = %d" % (len(peaks) + len(pits) - len(passes))
        print "len(crit_pts) = %d" % (len(peaks) + len(pits) + len(passes))
        if 1:
            ncontours = 30
            import pylab as pl
            # visualize(arr, crit_pts=cpts, cmap='gray', ncontours=ncontours, surf_network=contour_tree.to_undirected())
            # visualize(arr, crit_pts=cpts, cmap='gray', ncontours=ncontours, surf_network=None)
            # visualize(arr, cmap='gray', ncontours=ncontours)
            filled_arr = arr.copy()
            for region in sorted(regions, key=height_func, reverse=True):
                # filled_arr = arr.copy()
                X = [_[0] for _ in regions[region]]
                Y = [_[1] for _ in regions[region]]
                filled_arr[X,Y] = 2*arr.max()
                visualize(filled_arr, crit_pts=cpts, ncontours=None, cmap='gray', new_fig=False)
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

join_adj = {
        1: {},
        2: {1: {}},
        2.1: {2: {}},
        3: {2.1: {}},
        4: {3: {}},
        4.6: {4: {}},
        4.9: {4: {}},
        5: {4.6: {}},
        6: {4.9: {}},
        6.1: {5: {}},
        6.5: {5: {}},
        6.9: {6: {}},
        7: {6: {}},
        7.2: {6.1: {}},
        8: {6.9: {}},
        8.3: {7.2: {}},
        9: {6.5: {}},
        10: {8.3: {}},
        }

split_adj = {
        1: {2.1: {}},
        2: {3: {}},
        2.1: {3: {}},
        3: {4: {}},
        4: {4.6: {}},
        4.6: {4.9: {}},
        4.9: {5: {}},
        5: {6: {}},
        6: {6.1: {}},
        6.1: {6.5: {}},
        6.5: {6.9: {}},
        6.9: {7: {}},
        7: {7.2: {}},
        7.2: {8: {}},
        8: {8.3: {}},
        8.3: {9: {}},
        9: {10: {}},
        10: {}}
