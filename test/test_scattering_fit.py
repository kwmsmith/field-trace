from nose.tools import ok_, eq_, set_trace
import numpy as np
from field_trace.scattering_fit import dens_diff_perp, dens_diff_parallel

from field_trace.test.test_critical_point_network import random_periodic_upsample

class test_dens_diff(object):

    def setUp(self):
        self.nx, self.ny = 512, 512
        self.zz = np.zeros((self.nx, self.ny), dtype=np.double)

    def test_dd_perp_zeros(self):
        den_diff = dens_diff_perp(den=self.zz, idxs=range(10), outer_len=1.0)
        ok_(np.allclose(den_diff[:,1], 0.0))

    def test_dd_perp_const(self):
        den = self.zz
        den.fill(0.1)
        den_diff = dens_diff_perp(den=den, idxs=range(10), outer_len=1.0)
        ok_(np.allclose(den_diff[:,1], 0.0))

    def test_dd_perp_delta(self):
        den = self.zz
        den[self.nx/2, self.ny/2] = 1.0
        den_diff = dens_diff_perp(den=den, idxs=range(self.nx/8), outer_len=1.0)
        # print
        # print den_diff
        ok_(not np.allclose(den_diff[:,1], 0.0))

    def test_dd_parallel_zeros(self):
        den_diff = dens_diff_parallel(den=self.zz, idxs=[(0,0)], outer_len=1.0, nr=self.nx/4)
        ok_(np.allclose(den_diff[:,1], 0.0))

    def test_dd_parallel_const(self):
        den = self.zz
        den.fill(0.1)
        den_diff = dens_diff_parallel(den=den, idxs=[(0,0)], outer_len=1.0, nr=self.nx/4)
        ok_(np.allclose(den_diff[:,1], 0.0))

    def test_dd_parallel_delta(self):
        den = self.zz
        den[self.nx/2, self.ny/2] = 1.0
        den_diff = dens_diff_parallel(den=den, idxs=[(4,4)], outer_len=1.0, nr=self.nx/4)
        ok_(not np.allclose(den_diff[:,1], 0.0))

    def test_dd_parallel_random(self):
        arr = random_periodic_upsample(self.nx, 16, seed=3)
        nsamples = 300
        # idxs = [(0,0), (100,100), (70, 0), (256, 128)]
        idxs = zip(np.random.randint(0, self.nx, nsamples), np.random.randint(0, self.ny, nsamples))
        den_diff = dens_diff_parallel(den=arr, idxs=idxs, outer_len=self.nx, nr=self.nx/4)
        import pylab as pl
        pl.ion()
        pl.figure()
        pl.imshow(arr, interpolation='nearest', cmap='hot')
        pl.figure()
        pl.plot(den_diff[:,0], den_diff[:,1], 'ro-')
        pl.errorbar(den_diff[:,0], den_diff[:,1], yerr=den_diff[:,2])
        print
        print den_diff
        raw_input('enter to continue')
