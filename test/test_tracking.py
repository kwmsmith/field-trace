import tracking
import numpy as np

from nose.tools import ok_, eq_

class test_track(object):

    def setUp(self):
        self.npts = 500
        self.nx = self.ny = 100
        self.pts1 = np.random.randint(self.nx, size=(self.npts, 3))
        self.pts1[:,2] = np.random.randint(1, 10, size=self.npts)
        self.pts2 = np.random.randint(self.nx, size=(self.npts, 3))
        self.pts2[:,2] = np.random.randint(1, 10, size=self.npts)

    def _test_track_forward_backward(self):
        tracking.track_forward_backward(self.pts1, self.pts2, self.nx, self.ny)


    def test_track(self):
        map = tracking.track_regions(self.pts1, self.pts2, self.nx, self.ny)
        eq_(map, range(self.npts))
        self.pts2 = self.pts1.copy()
        self.pts2[0], self.pts2[1] = self.pts1[1], self.pts1[0]
        map = tracking.track_regions(self.pts1, self.pts2, self.nx, self.ny)
        check = range(self.npts)
        check[0], check[1] = check[1], check[0]
        eq_(map, check)
