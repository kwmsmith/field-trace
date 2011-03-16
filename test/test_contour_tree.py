from pprint import pprint, pformat
import networkx as nx

import contour_tree as ct

from nose.tools import eq_, ok_

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
        eq_(sorted(ct.peak_pit_nodes(join)), [7, 8, 9, 10])
        eq_(sorted(ct.peak_pit_nodes(split)), [1, 2])
        eq_(sorted(ct.pass_nodes(join)), [4, 5, 6])
        eq_(sorted(ct.pass_nodes(split)), [3])

        eq_(join.adj,
                {
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
                )

        eq_(split.adj,
                {
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
                    10: {}})

        if 1:
            import pylab as pl
            pl.ion()
            nx.draw_shell(join)
            pl.title('join')
            pl.figure()
            nx.draw_shell(split)
            pl.title('split')
            raw_input("enter to continue")

    def test_contour_tree(self):
        contour_tree = ct.contour_tree(self.mesh, self.height_func)
        import pylab as pl
        pl.ion()
        nx.draw(contour_tree)
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
