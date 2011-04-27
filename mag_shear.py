import _critical_points as _cp
from scipy.ndimage import gaussian_filter
import numpy as np
import pylab as pl
import tables
from kaw_analysis.vcalc import gradient
from region_analysis import radial_profiles

from test_tracking import h5fname

def mag_shear(psi_arr):
    r'''
    computes \grad_{\perp} \bvec{B} = 
             \frac{\grad{\psi}}{2 (\grad{\psi})^2} \cdot \grad{(\grad{\psi})^2}
    '''
    psi_grad_x, psi_grad_y = gradient(psi_arr)
    psi_grad_2 = psi_grad_x**2 + psi_grad_y**2
    psi_grad_grad_x, psi_grad_grad_y = gradient(psi_grad_2)
    mag_shear = psi_grad_x * psi_grad_grad_x + psi_grad_y * psi_grad_grad_y
    return mag_shear / (2.0 * psi_grad_2)

def vis_mag_shear(h5fname):
    dta = tables.openFile(h5fname, 'r')
    pl.figure()
    for idx, psi_arr in enumerate(dta.walkNodes('/psi', 'Array')):
        print psi_arr.name
        psi_arr = psi_arr.read()
        psi_arr = gaussian_filter(psi_arr, sigma=8.0, mode='wrap')
        mshear = mag_shear(psi_arr)
        pl.clf()
        pl.imshow(mshear, interpolation='nearest', cmap='hot')
        pl.savefig('mshear_%03d.pdf' % idx)
    dta.close()

def mag_shear_rad_profs(psi_arr):
    raise RuntimeError("do scatterplot rather than averaging...")
    psi_arr = gaussian_filter(psi_arr, sigma=2.0, mode='wrap')
    mshear = mag_shear(psi_arr)
    surf = _cp.TopoSurface(psi_arr)
    rprofs = radial_profiles(surf, threshold=25, other_arr=mshear)
    pl.ion()
    pl.figure()
    for minmax, (rprof, region) in rprofs.items():
        # minmax_flux = arr_div[minmax]
        pts, extremum_val, avg_shear, avg_shear_errs, avg_dists, avg_dists_errs = \
                zip(*rprof)
        # fluxes = np.abs(np.array(fluxes) - minmax_flux)
        # avg_fluxes = np.abs(np.array(avg_fluxes) - minmax_flux)
        # pl.plot(avg_dists, avg_fluxes, 'd-')
        pl.plot(avg_dists, avg_shear, 'd-')
        raw_input('enter to continue')
        pl.clf()
    import pdb; pdb.set_trace()
    pl.close('all')

def save_mshear_rad_profs(h5fname):
    dta = tables.openFile(h5fname, 'r')
    for idx, psi_arr in enumerate(dta.walkNodes('/psi', 'Array')):
        print psi_arr.name
        psi_arr = psi_arr.read()
        mag_shear_rad_profs(psi_arr)
    dta.close()

save_mshear_rad_profs('./data.h5')
