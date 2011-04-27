from scipy import stats
import _critical_points as _cp
from region_analysis import radial_profiles

from kaw_analysis import vcalc

import numpy as np

# from nose.tools import ok_, eq_

from test_critical_point_network import random_periodic_upsample

import pylab as pl
pl.ion()

def test_radial_profiles():
    arr = random_periodic_upsample(128, 16, seed=0)
    mask = np.zeros(arr.shape, dtype=np.bool_)
    arr_x = vcalc.cderivative(arr, 'X_DIR')
    arr_y = vcalc.cderivative(arr, 'Y_DIR')
    arr_div = np.sqrt(arr_x**2 + arr_y**2)
    surf = _cp.TopoSurface(arr)
    rprofs = radial_profiles(surf, threshold=25, expand_regions=1, other_arr=arr_div, mask=mask)
    arr[mask] = 2 * arr.max()
    pl.imshow(arr, interpolation='nearest')
    pl.figure()
    pl.imshow(arr_div)
    pl.figure()
    pl.hold(True)
    linreg_xy = ([], [])
    for minmax, (rprof, region) in rprofs.items():
        # minmax_flux = arr_div[minmax]
        pts, fluxes, avg_fluxes, avg_fluxes_errs, avg_dists, avg_dists_errs = \
                zip(*rprof)
        linreg_xy[0].extend(fluxes)
        linreg_xy[1].extend(avg_fluxes)
        # fluxes = np.abs(np.array(fluxes) - minmax_flux)
        # avg_fluxes = np.abs(np.array(avg_fluxes) - minmax_flux)
        # pl.plot(avg_dists, avg_fluxes, 'd-')
        pl.plot(avg_dists, avg_fluxes, 'd-')
    pl.grid()
    slope, intercept, rval, pval, stderr = stats.linregress(*linreg_xy)
    print
    print "slope: %f" % slope
    print "intercept: %f" % intercept
    print "rval: %f" % rval
    print "pval: %f" % pval
    print "stderr: %f" % stderr
    import pdb; pdb.set_trace()
