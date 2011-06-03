import numpy as np
from map_reduce import map_reduce
from mag_shear import integrate_theta

def dens_diff_perp(den, idxs, outer_len):
    '''
    integrates perpendicular to the 2D density field.
    '''
    nx, ny = den.shape
    assert nx == ny
    dz = outer_len / nx
    res = []
    for idx in idxs:
        x0_pencil = den[idx, :]
        for j in range(nx):
            delta_x = np.abs(idx - j)
            x1_pencil = den[j, :]
            ss = np.sum(x0_pencil - x1_pencil) * dz
            res.append((delta_x % (nx/2), ss))
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
