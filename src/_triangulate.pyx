cdef extern from "math.h":
    double sqrt(double)

cdef class vec:
    cdef double x, y, z

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    property x:
        def __get__(self):
            return self.x

    property y:
        def __get__(self):
            return self.y

    property z:
        def __get__(self):
            return self.z

    cpdef vec crossp(s, o):
        return vec(s.y * o.z - o.y * s.z,
                   -s.x * o.z + o.x * s.z,
                   s.x * o.y - o.x * s.y)

    def __richcmp__(self, other, int op):
        if op == 2: # eqiv to __eq__(self, other)
            return (self.x == other.x and
                    self.y == other.y and
                    self.z == other.z)
        else:
            raise NotImplementedError()

    cpdef vec sub(self, other):
        return vec(self.x - other.x,
                   self.y - other.y,
                   self.z - other.z)

    def __sub__(self, other):
        return self.sub(other)

    def __rsub__(self, other):
        return other.sub(self)

    cpdef double dot(self, other):
        return (self.x * other.x +
                self.y * other.y +
                self.z * other.z)

    def norm(self):
        return sqrt(self.dot(self))

def tri_plane_cos(vec a, vec b, vec c, vec d):
    cdef vec v1 = b.sub(a)
    cdef vec v2 = d.sub(a)
    cdef vec v3 = b.sub(c)
    cdef vec v4 = d.sub(c)
    cdef vec n1 = v2.crossp(v1)
    cdef vec n2 = v3.crossp(v4)
    return n1.dot(n2) / (n1.norm() * n2.norm())

def tri_plane_cos_from_z(za, zb, zc, zd):
    # A ----- B
    # |       |
    # |       |
    # |       |
    # D ----- C
    return tri_plane_cos(vec(0, 0, z=za),
                         vec(0, 1, z=zb),
                         vec(1, 1, z=zc),
                         vec(1, 0, z=zd))

class graph:

    def __init__(self):
        from collections import defaultdict
        self._g = defaultdict(set)

    def add_edge(self, a, b):
        if b not in self._g[a]:
            self._g[a].add(b)
            self._g[b].add(a)

def mesh(arr):
    G = graph()
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
            tpc_bd = tri_plane_cos_from_z(a, b, c, d)
            tpc_ac = tri_plane_cos_from_z(b, c, d, a)
            if tpc_bd >= tpc_ac:
                G.add_edge((i,(j+1)%ny), ((i+1)%nx, j))
            else:
                G.add_edge((i,j), ((i+1)%nx, (j+1)%ny))
    return G
