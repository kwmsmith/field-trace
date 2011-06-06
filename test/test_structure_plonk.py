import numpy as np
from numpy import random
from field_trace.structure_plonk import structure_plonk, gaus_fit, gaussian

import pylab as pl

def _test_structure_plonk():
    arr = structure_plonk(size=512, nstruts=40, amp=2.0, rad_exp=3, core_radius=5, seed=1)
    pl.ion()
    pl.imshow(arr, interpolation='nearest', cmap='hot')
    raw_input('enter to continue')

def test_leastsq():
    pl.ion()
    dta = random.normal(loc=2.0, scale=3.0, size=10000)
    dta_hist = np.histogram(dta, normed=True, bins=np.sqrt(dta.size))
    x0, sig0 = dta.mean(), dta.std()
    best_x, best_sig = gaus_fit(dta_hist, x0, sig0)
    pl.figure()
    n, bins, patches = pl.hist(dta, normed=True, bins=np.sqrt(dta.size), histtype='step')
    fit = gaussian(bins, best_x, best_sig)
    pl.plot(bins, fit)
    raw_input('enter to continue')
