from wrap_gsl_interp2d import Interp2DPeriodic
from wrap_trace_integrator import TraceIntegrator

import numpy as np

def test_ti():

    # number of gridpoints.
    nn = 256
    # wavenumbers in X & Y dimensions, resp.
    n, m = 10, 3 
    x0max = x1max = 2.0 * np.pi
    Y = np.linspace(0, x1max, nn, endpoint=False)
    # print Y
    X = Y[:,np.newaxis]
    arr = np.sin(2. * np.pi * n * X / x0max) + np.sin(2. * np.pi * m * Y / x1max)

    interp = Interp2DPeriodic(x0max, x1max, arr)

    ntracers = 10
    scale = 1.0
    ti = TraceIntegrator(ntracers, scale, x0max, x1max, interp, interp, 1.e-6, 0.0)

    print "here test_ti"
