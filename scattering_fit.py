import numpy as np
from map_reduce import map_reduce
from mag_shear import integrate_theta

from contour_tree import wraparound_dist_1d

import pylab as pl

from structure_plonk import structure_plonk

def dens_diff_perp(den, idxs, outer_len):
    '''
    integrates perpendicular to the 2D density field.
    '''
    outer_len = float(outer_len)
    nx, ny = den.shape
    wdist = wraparound_dist_1d(nx)
    # assert nx == ny
    dz = outer_len / ny
    res = []
    for idx in idxs:
        x0_pencil = den[idx, :]
        for j in range(nx):
            delta_x = wdist(idx, j)
            x1_pencil = den[j, :]
            ss = np.sum(x0_pencil - x1_pencil) * dz
            res.append((delta_x, ss))
    def mapper(elm):
        return elm[0], elm[1]
    def reducer(gp):
        return sum(gp) / len(gp)
    dens_diff = map_reduce(res, mapper, reducer)
    dd_arr = np.array(sorted(dens_diff.items()))
    dd_arr[:,0] *= dz
    return dd_arr

def dens_diff_parallel(den, idxs, outer_len, nr=None):
    nx, ny = den.shape
    nr = nr or nx
    assert nx == ny
    dz = outer_len / nx
    center_x, center_y = nx/2, ny/2
    reduced_diffs = []
    for idx in idxs:
        ix, iy = idx
        den0 = den[idx]
        delta_den = -den + den0
        res = integrate_theta(delta_den, ix, iy, nr=nr)
        reduced_diffs.append(res)
    rs = reduced_diffs[0]['r'] * dz
    reduced_vals = np.array([rd['val'] for rd in reduced_diffs])
    vals = np.mean(reduced_vals, axis=0)
    valstds = np.std(reduced_vals, axis=0)
    dd_arr = np.empty((len(rs), 3), dtype=np.double)
    dd_arr[:,0] = rs
    dd_arr[:,1] = vals
    dd_arr[:,2] = valstds / np.sqrt(len(reduced_vals))
    return dd_arr

def scattering_fit_sim_structures(size, repeat, nstruts, rad_exp, seed=None):
    np.random.seed(seed)
    arrs = [structure_plonk(size=size, nstruts=nstruts, amp=1.0, rad_exp=rad_exp, core_radius=5) \
            for _ in range(repeat)]
    arr = np.concatenate(arrs, axis=1)
    idxs = range(0, size, max(size/16,1))
    pl.ion()
    pl.figure()
    for idx in idxs:
        print idx
        dd_arr = dens_diff_perp(arr, [idx], outer_len=repeat)
        pl.plot(dd_arr[:,0], dd_arr[:,1], 'o-')
    raw_input('enter to continue')

if __name__ == '__main__':
    scattering_fit_sim_structures(size=512, repeat=10, nstruts=512, rad_exp=1, seed=1)
