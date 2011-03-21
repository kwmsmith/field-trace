from collections import defaultdict

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

class Graph(object):

    def __init__(self):
        self._adj = defaultdict(set)

    def neighbors(self, node):
        if node in self._adj:
            return self._adj[node]
        else:
            raise KeyError(str(node))

    neighbors_iter = neighbors

    def degree(self, node):
        if node in self._adj:
            return len(self._adj[node])
        else:
            raise KeyError(str(node))

    def nodes(self):
        return self._adj.keys()

    def add_edge(self, a, b):
        self._adj[a].add(b)
        self._adj[b].add(a)
        assert b in self._adj[a]
        assert a in self._adj[b]

    def add_edges_from(self, edges):
        for (a, b) in edges:
            self.add_edge(a, b)

    def remove_edge(self, a, b):
        if a in self._adj and b in self._adj:
            self._adj[a].remove(b)
            if a != b:
                self._adj[b].remove(a)
        else:
            raise KeyError("edge (%s, %s) not in graph." % (a, b))

    def remove_node(self, node):
        if node in self._adj:
            for nbr in self.neighbors(node).copy():
                self.remove_edge(node, nbr)
            del self._adj[node]
        else:
            raise KeyError(str(node))

    def __len__(self):
        return len(self._adj)

    def edges(self):
        edges = set()
        for start in self._adj:
            for next in self.neighbors(start):
                if (next, start) not in edges:
                    edges.add((start, next))
        return edges

class DiGraph(object):

    def __init__(self):
        self._adj = defaultdict(set)
        self._inv = defaultdict(set)

    def successors(self, node):
        if node in self._adj:
            return self._adj[node]
        else:
            raise KeyError(str(node))

    def predecessors(self, node):
        if node in self._inv:
            return self._inv[node]
        else:
            raise KeyError(str(node))

    def in_degree(self, node=None):
        if node is None:
            return dict([(n, self._in_degree(n)) for n in self._adj])
        else:
            return self._in_degree(node)

    def _in_degree(self, node):
        if node in self._inv:
            return len(self._inv[node])
        else:
            raise KeyError(str(node))

    def out_degree(self, node=None):
        if node is None:
            return dict([(n, self._out_degree(n)) for n in self._adj])
        else:
            return self._out_degree(node)

    def _out_degree(self, node):
        if node in self._adj:
            return len(self._adj[node])
        else:
            raise KeyError(str(node))

    def add_edge(self, a, b):
        self._adj[a].add(b)
        self._adj[b]
        self._inv[b].add(a)
        self._inv[a]

    def add_edges_from(self, edges):
        for (a, b) in edges:
            self.add_edge(a, b)

    def remove_edge(self, a, b):
        if a in self._adj and b in self._inv:
            self._adj[a].remove(b)
            self._inv[b].remove(a)
        else:
            raise KeyError("%s or %s not in graph." % (a, b))

    def remove_node(self, node):
        if node in self._adj:
            for succ in self.successors(node).copy():
                self.remove_edge(node, succ)
            for pred in self.predecessors(node).copy():
                self.remove_edge(pred, node)
            del self._adj[node]
            if node in self._inv:
                del self._inv[node]
        else:
            raise KeyError(str(node))

    def __len__(self):
        return len(self._adj) + len(self._inv)

    order = __len__

    def nodes(self):
        return self._adj.keys()

    def edges(self):
        edges = []
        for start in self._adj:
            for next in self.successors(start):
                edges.append((start, next))
        return edges
