import numpy as np
from numpy.fft import rfft2, irfft2
from kaw_analysis import test_vcalc as tv
from kaw_analysis import vcalc

from wrap_gsl_interp2d import Interp2DPeriodic

from nose.tools import ok_, eq_, set_trace

import field_trace

global_N = 64

import pylab as pl

import upsample

def test_upsample():
    upsample_factor = 4
    N = global_N
    test_data = tv.sin_cos_prod(N, N/4, N/8)
    test_data_4 = tv.sin_cos_prod(N*4, N/4, N/8)

    upsampled = upsample.upsample(test_data, factor=4)

    ok_(np.allclose(upsample_factor**2 * upsampled, test_data_4))

    if 0:
        pl.ion()
        pl.imshow(test_data, interpolation='nearest', cmap='hot')
        pl.title('original data')
        pl.figure()
        pl.imshow(upsampled, interpolation='nearest', cmap='hot')
        pl.title('upsampled by %d' % upsample_factor)
        pl.figure()
        pl.imshow(test_data_4, interpolation='nearest', cmap='hot')
        pl.title('comparison')
        raw_input('enter to continue')
