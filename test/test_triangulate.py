import numpy as np
import _triangulate as _tri

from kaw_analysis import test_vcalc as tv

from nose.tools import eq_, ok_, set_trace

def test_vec():
    v = _tri.vec(1,2,3)
    eq_(v.x, 1.0)
    eq_(v.y, 2.0)
    eq_(v.z, 3.0)

def test_crossp():
    # trivial case
    eq_(_tri.vec(0,0,0), _tri.vec().crossp(_tri.vec(10,20,30)))
    # test parallel vectors
    eq_(_tri.vec(0,0,0), _tri.vec(1,2,3).crossp(_tri.vec(10,20,30)))
    # test anti-parallel vectors
    eq_(_tri.vec(0,0,0), _tri.vec(-1,-2,-3).crossp(_tri.vec(10,20,30)))
    # test perp vectors
    eq_(_tri.vec(z=1), _tri.vec(x=1).crossp(_tri.vec(y=1)))
    # test perp vectors
    eq_(_tri.vec(z=-1), _tri.vec(y=1).crossp(_tri.vec(x=1)))

def test_sub():
    v1 = _tri.vec(1,2,3)
    v2 = _tri.vec(10, 13, 21)
    eq_(_tri.vec(9, 11, 18), v2.__sub__(v1))
    eq_(_tri.vec(9, 11, 18), v2 - v1)
    eq_(_tri.vec(-9, -11, -18), v1 - v2)
    eq_(_tri.vec(-9, -11, -18), v1.__sub__(v2))

def test_dotp():
    v1 = _tri.vec()
    v2 = _tri.vec(1,2,3)
    eq_(0.0, v1.dot(v2))
    eq_(1 + 2**2 + 3**2, v2.dot(v2))
    eq_(0.0, v1.norm())
    eq_(5.0, _tri.vec(y=5).norm())

def test_tri_plane_cos():
    a = 1.0
    b = 0.0
    c = 1.0
    d = 0.0
    tpc1 = _tri.tri_plane_cos_from_z(a, b, c, d)
    tpc2 = _tri.tri_plane_cos_from_z(b, c, d, a)
    eq_(tpc1, tpc2)
    eq_(1.0, _tri.tri_plane_cos_from_z(*[0.0]*4))
    a = b = d = 0.0
    lim = 10000
    calc = 1. / (np.sqrt(2) * np.sqrt(0.5 + np.arange(lim)**2))
    tpcs = [_tri.tri_plane_cos_from_z(a, b, c, d) for c in range(lim)]
    tpcs2 = [_tri.tri_plane_cos_from_z(b, c, d, a) for c in range(lim)]
    ok_(np.allclose(tpcs, calc))

def test_triangulate():
    arr = tv.sin_cos_prod(20, 8, 3)
    gr = _tri.mesh(arr)
    visualize_mesh(arr, gr)
    raw_input('enter to continue')

def visualize_mesh(arr, gr):
    import pylab as pl
    pl.ion()
    # pl.imshow(arr, interpolation='nearest', cmap='jet')
    pl.imshow(arr, cmap='jet')
    nx, ny = arr.shape
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            other_pts = gr._g[i,j]
            for x,y in other_pts:
                pl.plot([j,y], [i,x], 'k--')
