from wrap_gsl_interp2d import Interp2DPeriodic
from wrap_trace_integrator import TraceIntegrator

import numpy as np
import pylab as pl

from kaw_analysis import vcalc

def test_ti():

    # number of gridpoints.
    nn = 256
    # wavenumbers in X & Y dimensions, resp.
    n, m = 20, 13 
    x1max = x0max =  256.0
    Y = np.linspace(0, x1max, nn, endpoint=False)
    X = Y[:,np.newaxis]
    arr = np.sin(2. * np.pi * n * X / x0max) + np.sin(2. * np.pi * m * Y / x1max)

    d_x = vcalc.cderivative(arr, 'X_DIR')
    d_y = vcalc.cderivative(arr, 'Y_DIR')
    deriv_norm = 10.0 * 1.0 / np.sqrt((d_x**2 + d_y**2).max())

    d_x *= deriv_norm
    d_y *= deriv_norm

    interp1 = Interp2DPeriodic(x0max, x0max, -d_y)
    interp2 = Interp2DPeriodic(x0max, x0max, d_x)

    ntracers = 10
    scale = 1.0
    ti = TraceIntegrator(ntracers, scale, x0max, x1max, interp1, interp2, 1.e-6, 0.0)

    positions = np.empty(2*ntracers, dtype=np.double)
    positions[::2] = np.linspace(0.0, nn, ntracers)
    positions[1::2] = np.linspace(0.0, nn, ntracers)

    times = np.linspace(0.0, 50.0, 10)

    positions_vs_t = []

    t_init = times[0]
    for time in times[1:]:
        positions_vs_t.append(list(positions.copy()))
        ti.evolve(t=t_init, t1=time, h=1.0e-6, y=positions)
        t_init = time

    if 0:
        pl.ion()
        dta = np.array(positions_vs_t)
        pl.imshow(arr)
        for i in range(0, 2*ntracers, 2):
            x0 = dta[:,i]
            y0 = dta[:,i+1]
            pl.scatter(y0, x0)
        raw_input("enter to quit")
