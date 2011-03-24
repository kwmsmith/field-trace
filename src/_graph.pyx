# cython: profile=False

# successors -- check
# predecessors -- check
# neighbors -- check
# in_degree -- check
# out_degree -- check
# degree -- check

# add_edge
# add_edges_from
# remove_edge

# nodes

cdef class Graph(object):

    cdef readonly dict adj

    def __init__(self):
        self.adj = {}

    def neighbors(self, node):
        # try:
        return list(self.adj[node])
        # except KeyError:
            # raise KeyError("node %s is not in the graph" % (node,))

    def neighbors_iter(self, node):
        # try:
        return list(self.adj[node])
        # except KeyError:
            # raise KeyError("node %s is not in the graph" % (node,))

    def degree_iter(self, node=None):
        if node is None:
            nodes_nbrs = iter(self.adj.items())
        else:
            nodes_nbrs=[(n,self.adj[n]) for n in self.adj[node]]
        return [(n, len(nbrs)+(n in nbrs)) for n,nbrs in nodes_nbrs]

    def degree(self, node=None):
        if node in self.adj: # return a single node
            return self.degree_iter(node)[1]
        else:           # return a dict
            return dict(self.degree_iter(node))

    def nodes_iter(self):
        return iter(self.adj.keys())

    def nodes(self):
        return list(self.nodes_iter())

    def add_edge(self, u, v):
        # add nodes
        if u not in self.adj:
            self.adj[u] = set()
        if v not in self.adj:
            self.adj[v] = set()
        self.adj[u].add(v)
        self.adj[v].add(u)

    def add_edges_from(self, edges):
        for (a, b) in edges:
            self.add_edge(a, b)

    def remove_edge(self, u, v):
        # try:
        self.adj[u].remove(v)
        if u != v:  # self-loop needs only one entry removed
            self.adj[v].remove(u)
        # except KeyError:
            # raise KeyError("The edge %s-%s is not in the graph"%(u,v))

    def remove_node(self, n):
        adj = self.adj
        # try:
        nbrs = adj[n] # keys handles self-loops (allow mutation later)
        # except KeyError: # NetworkXError if n not in self
            # raise KeyError("The node %s is not in the graph."%(n,))
        if n not in nbrs:
            for u in nbrs:
                adj[u].remove(n)   # remove all edges n-u in graph
        else:
            for u in nbrs.copy():
                adj[u].remove(n)   # remove all edges n-u in graph
        del adj[n]          # now remove node

    def __len__(self):
        return len(self.adj)

    order = __len__

    def edges_iter(self):
        seen=set()     # helper dict to keep track of multiply stored edges
        edges = []
        nodes_nbrs = iter(self.adj.items())
        for n,nbrs in nodes_nbrs:
            for nbr in nbrs:
                if nbr not in seen:
                    edges.append((n,nbr))
            seen.add(n)
        return edges
    
    def edges(self):
        return self.edges_iter()

cdef class DiGraph(Graph):

    cdef readonly dict pred, succ

    def __init__(self):
        self.adj = {}  # empty adjacency dictionary
        self.pred = {}  # predecessor
        self.succ = self.adj  # successor

    def successors(self, n):
        # try:
        return list(self.succ[n])
        # except KeyError:
            # raise KeyError("The node %s is not in the digraph."%(n,))

    def predecessors(self, n):
        # try:
        return list(self.pred[n])
        # except KeyError:
            # raise KeyError("The node %s is not in the digraph."%(n,))

    def in_degree_iter(self, node=None):
        if node is None:
            nodes_nbrs = self.pred.items()
        else:
            nodes_nbrs=[(node,self.pred[node])]
        return [(n, len(nbrs)) for n,nbrs in nodes_nbrs]

    def in_degree(self, nbunch=None):
        if nbunch in self.adj:      # return a single node
            return self.in_degree_iter(nbunch)[0][1]
        else:           # return a dict
            return dict(self.in_degree_iter(nbunch))

    def out_degree_iter(self, node=None):
        if node is None:
            nodes_nbrs= self.succ.items()
        else:
            nodes_nbrs=[(node,self.succ[node])]
        return [(n, len(nbrs)) for n,nbrs in nodes_nbrs]

    def out_degree(self, nbunch=None):
        if nbunch in self.adj:      # return a single node
            return self.out_degree_iter(nbunch)[0][1]
        else:           # return a dict
            return dict(self.out_degree_iter(nbunch))

    def add_edge(self, u, v):
        # add nodes
        if u not in self.succ:
            self.succ[u] = set()
            self.pred[u] = set()
        if v not in self.succ:
            self.succ[v] = set()
            self.pred[v] = set()
        # add the edge
        self.succ[u].add(v)
        self.pred[v].add(u)

    def add_edges_from(self, edges):
        for (a, b) in edges:
            self.add_edge(a, b)

    def remove_edge(self, u, v):
        # try:
        self.succ[u].remove(v)
        self.pred[v].remove(u)
        # except KeyError:
            # raise KeyError("The edge %s-%s not in graph."%(u,v))

    def remove_node(self, n):
        # try:
        nbrs=self.succ[n]
        # except KeyError: # NetworkXError if n not in self
            # raise KeyError("The node %s is not in the digraph."%(n,))
        for u in nbrs:
            self.pred[u].remove(n)
        del self.succ[n]          # remove node from succ
        for u in self.pred[n]:
            self.succ[u].remove(n)
        del self.pred[n]          # remove node from pred

    def edges_iter(self):
        return [(n, nbr) for n,nbrs in self.adj.items() for nbr in nbrs]
