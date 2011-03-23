import networkx as nx

import _graph as gr

from random import randrange, shuffle

from nose.tools import ok_, eq_, set_trace

N = 50

def compare_graphs(g1, g2):
    eq_(set(g1.nodes()), set(g2.nodes()))
    eq_(set(g1.edges()), set(g2.edges()))

def test_random():
    nx_g = nx.Graph()
    gr_g = gr.Graph()
    nx_dg = nx.DiGraph()
    gr_dg = gr.DiGraph()
    for _ in range(N):
        for i in range(N):
            a = randrange(N)
            b = randrange(N)
            nx_dg.add_edge(a, b)
            nx_g.add_edge(a, b)
            gr_dg.add_edge(a, b)
            gr_g.add_edge(a, b)
        compare_graphs(gr_dg, nx_dg)
        compare_graphs(gr_g, nx_g)

def test_random_removes():
    nx_g = nx.Graph()
    gr_g = gr.Graph()
    nx_dg = nx.DiGraph()
    gr_dg = gr.DiGraph()
    for _ in range(N):
        for i in range(N**2):
            a = randrange(N)
            b = randrange(N)
            nx_dg.add_edge(a, b)
            nx_g.add_edge(a, b)
            gr_dg.add_edge(a, b)
            gr_g.add_edge(a, b)
        for e in gr_dg.edges():
            gr_dg.remove_edge(*e)
            nx_dg.remove_edge(*e)
        compare_graphs(gr_dg, nx_dg)
        ok_(not gr_dg.edges())
        ok_(not nx_dg.edges())
        for e in gr_g.edges():
            gr_g.remove_edge(*e)
            nx_g.remove_edge(*e)
        ok_(not gr_g.edges())
        ok_(not nx_g.edges())
        compare_graphs(gr_g, nx_g)

def test_random_remove_nodes():
    nx_g = nx.Graph()
    gr_g = gr.Graph()
    nx_dg = nx.DiGraph()
    gr_dg = gr.DiGraph()
    for _ in range(N):
        for i in range(N**2):
            a = randrange(N)
            b = randrange(N)
            nx_dg.add_edge(a, b)
            nx_g.add_edge(a, b)
            gr_dg.add_edge(a, b)
            gr_g.add_edge(a, b)
        for n in gr_dg.nodes():
            gr_dg.remove_node(n)
            nx_dg.remove_node(n)
        compare_graphs(gr_dg, nx_dg)
        ok_(not gr_dg.nodes())
        ok_(not nx_dg.nodes())
        ok_(not gr_dg.edges())
        ok_(not nx_dg.edges())
        for n in gr_g.nodes():
            gr_g.remove_node(n)
            nx_g.remove_node(n)
        ok_(not gr_g.edges())
        ok_(not nx_g.edges())
        ok_(not gr_g.nodes())
        ok_(not nx_g.nodes())
        compare_graphs(gr_g, nx_g)
