from itertools import izip
import field_trace
from numpy import pi
import _critical_points as _cp
from scipy.ndimage import gaussian_filter
import numpy as np
import pylab as pl
import tables
from kaw_analysis.vcalc import gradient
from region_analysis import radial_profiles, flux_tube_radial_scatter, expand_region_circ, flux_tube_radial_spokes
from contour_tree import wraparound_dist_vec

# from test_tracking import h5fname

SIGMA = 2.0

def total_shear(bx, by, x0, y0):
    raise RuntimeError("doesn't work; need better treatment of gradients near "
            "R=0.  Currently get wild oscillations near R=0, since gradient uses "
            "fourier transforms.")
    nx, ny = bx.shape
    br, btheta = cartesian_to_polar(bx, by, x0, y0)
    X, Y, R, theta = cartesian_to_polar_coords(x0, y0, nx, ny)
    br_by_r = br/R
    br_by_r[x0, y0] = 0.0
    br_by_r_x, br_by_r_y = gradient(br_by_r)
    bt_by_r = btheta/R
    bt_by_r[x0, y0] = 0.0
    btheta_by_r_x, btheta_by_r_y = gradient(bt_by_r)
    shear_1, _ = cartesian_to_polar(btheta_by_r_x, btheta_by_r_y, x0, y0)
    _, shear_2 = cartesian_to_polar(br_by_r_x, br_by_r_y, x0, y0)
    return shear_1 + shear_2

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
        psi_arr = gaussian_filter(psi_arr, sigma=SIGMA, mode='wrap')
        mshear = mag_shear(psi_arr)
        pl.clf()
        pl.imshow(mshear, interpolation='nearest', cmap='hot')
        pl.savefig('mshear_%03d.pdf' % idx)
    dta.close()

def rad_scatter_driver(h5fname):
    dta = tables.openFile(h5fname, 'r')
    walkers = [dta.walkNodes('/psi', 'Array'),
               dta.walkNodes('/cur', 'Array'),
               dta.walkNodes('/den', 'Array')]
    for idx, arrs in enumerate(izip(*walkers)):
        if idx < 100: continue
        psi_arr, cur_arr, den_arr = arrs
        print psi_arr.name
        psi_arr = psi_arr.read()
        psi_arr = psi_arr.astype(np.double)
        den_arr = den_arr.read()
        den_arr = den_arr.astype(np.double)
        cur_arr = cur_arr.read()
        cur_arr = cur_arr.astype(np.double)
        # psi_arr = gaussian_filter(psi_arr, sigma=SIGMA, mode='wrap')
        # mag_shear_rad_scatter(psi_arr)
        mag_shear_theta_sectors(psi_arr, cur_arr, den_arr)
        # mag_shear_rad_spokes(psi_arr)
    dta.close()

def cartesian_to_polar_coords(x0, y0, nx, ny):
    X, Y = np.ogrid[0:nx, 0:ny]
    X = X - x0
    Y = Y - y0
    xsgn = np.sign(X)
    ysgn = np.sign(Y)
    X = np.where(np.abs(X) > nx/2, X - xsgn * nx, X)
    Y = np.where(np.abs(Y) > ny/2, Y - ysgn * ny, Y)
    R = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)
    theta = (theta + 2*pi) % (2*pi)
    return X, Y, R, theta

def cartesian_to_polar(ax, ay, x0, y0):
    nx, ny = ax.shape
    X, Y, R, theta = cartesian_to_polar_coords(x0, y0, nx, ny)
    ar = (ax * X + ay * Y) / R
    atheta = (- ax * Y + ay * X) / R
    atheta[x0, y0] = 0.0
    return (ar, atheta)

def mag_shear_rad_spokes(psi_arr):
    nx, ny = psi_arr.shape
    wdist = wraparound_dist_vec(nx, ny)
    # mshear = mag_shear(psi_arr)
    psi_grad_x, psi_grad_y = gradient(psi_arr)
    bmag = np.sqrt(psi_grad_x**2 + psi_grad_y**2)
    surf = _cp.TopoSurface(psi_arr)
    regions = surf.get_minmax_regions()
    pl.ion()
    pl.figure()
    for (minmax, pss, type), region in regions.items():
        br, btheta = cartesian_to_polar(-psi_grad_y, psi_grad_x, minmax[0], minmax[1])
        if len(region) < 100:
            continue
        region = expand_region_circ(psi_arr, region, minmax, extra=10.0)
        spokes = flux_tube_radial_spokes(region, minmax, psi_arr, type)
        pl.clf()
        # max_2_bmags = []
        max_2_btheta = []
        for spoke in spokes:
            sp_arr = np.array(spoke)
            xs, ys = sp_arr[:,0], sp_arr[:,1]
            # bmags = bmag[xs, ys]
            bthetas = btheta[xs, ys]
            max_2_btheta.append((bthetas.max(), bthetas, spoke))
        max_2_btheta.sort(key=lambda x: x[0], reverse=True)
        spokes = [sp for (m, bthetas, sp) in max_2_btheta]
        btheta_cpy = np.zeros_like(btheta)
        bmag_cpy = np.zeros_like(bmag)
        br_cpy = np.zeros_like(br)
        xs, ys = zip(*region)
        btheta_cpy[xs, ys] = btheta[xs, ys]
        bmag_cpy[xs, ys] = bmag[xs, ys]
        br_cpy[xs, ys] = br[xs, ys]
        for spoke in spokes:
            sp_arr = np.array(spoke)
            xs, ys = sp_arr[:,0], sp_arr[:,1]
            dists = wdist(minmax[0], minmax[1], xs, ys)
            bthetas = btheta[xs, ys]
            # mshears = mshear[xs, ys]
            pl.subplot(221)
            pl.plot(dists, bthetas, 'o-')
            # pl.subplot(222)
            # pl.plot(dists, mshears, 's-')
            pl.subplot(222)
            pl.plot(dists[1:], bthetas[1:]/dists[1:], 'd-')
        pl.subplot(221)
        pl.grid()
        pl.title(r'$B_{\theta}$ vs. $r$')
        # pl.subplot(222)
        # pl.grid()
        # pl.title(r'$\nabla_{\perp} B$ vs. $r$')
        pl.subplot(222)
        pl.grid()
        pl.title(r'$B_{\theta}/r$ vs. $r$')
        pl.subplot(234)
        pl.imshow(btheta_cpy, interpolation='nearest', cmap='hot')
        pl.title(r'$B_{\theta}$')
        pl.subplot(235)
        pl.imshow(br_cpy, interpolation='nearest', cmap='hot')
        pl.title(r'$B_{r}$')
        pl.subplot(236)
        pl.imshow(bmag_cpy, interpolation='nearest', cmap='hot')
        pl.title(r'$|B|$')
        raw_input('enter to continue')
        pl.clf()
    pl.close('all')

def theta_sectors(arr, x0, y0, region, nr=None, ntheta=None):
    nx, ny = arr.shape
    X, Y, R, theta = cartesian_to_polar_coords(x0, y0, nx, ny)
    idx0, idx1 = zip(*region)
    # idx0, idx1 = np.where(R<=rmax)
    Rs = R[idx0, idx1]
    rmax = Rs.max()
    thetas = theta[idx0, idx1]
    vals = arr[idx0, idx1]
    nr = nr or int(rmax+1)
    ntheta = ntheta or 1
    rbins = (Rs / rmax * (nr - 1)).astype(np.int32)
    thetabins = (thetas / (2*pi) * (ntheta -1)).astype(np.int32)
    dta = zip(rbins, thetabins, Rs, thetas, vals)
    def mapper(elm):
        rbin, thetabin, r, theta, v = elm
        return (rbin, thetabin), elm
    def reducer(gp):
        dta = np.array(gp,
                dtype=[('rbin', np.int32), ('thetabin', np.int32), ('r', np.float), ('theta', np.float), ('val', np.float)])
        rs = dta['r']
        thetas = dta['theta']
        vals = dta['val']
        return (rs.mean(), rs.std(), thetas.mean(), thetas.std(), vals.mean(), vals.std())
    from map_reduce import map_reduce
    sectors = map_reduce(dta, mapper, reducer)
    theta_sectors = {}
    for rbin, tbin in sectors:
        theta_sectors.setdefault(tbin, []).append(sectors[(rbin, tbin)])
    for tbin in theta_sectors:
        gp = theta_sectors[tbin]
        dta = np.array(gp,
                dtype=[('r', np.float), ('rstd', np.float),
                       ('theta', np.float), ('thetastd', np.float),
                       ('val', np.float), ('valstd', np.float)])
        dta.sort(order=['r'])
        theta_sectors[tbin] = dta
    return theta_sectors.values()

# for these fields, do an image plot & a radial profile plot:

# \psi
# J
# B_theta
# B_shear
# Gaussian Curvature
# n
# grad n
# |V|

def plot_theta_sectors(ax, field, minmax, region, ntheta, title, xvis=False):
    th_sectors = theta_sectors(field, minmax[0], minmax[1], region, ntheta=ntheta)
    ax.axhline(y=0, color='orange', linewidth=4)
    for sect in th_sectors:
        ax.plot(sect['r'], sect['val'], 's-')
        ax.errorbar(sect['r'], sect['val'],
                xerr=sect['rstd'], yerr=sect['valstd'])
    ax.grid()
    if not xvis:
        ax.tick_params(labelbottom=False)
    # ax.title(title)
    return th_sectors

def plot_sect_derivs_by_r(sectors, title):
    for sect in sectors:
        delta_r = np.diff(sect['r'])
        delta_btheta_by_r = np.diff(sect['val']/sect['r'])
        deriv = delta_btheta_by_r / delta_r
        pl.plot(sect['r'][:-1]+delta_r/2, deriv, 'o-')
    pl.grid()
    pl.title(title)

def field_region_image(ax, field, region, minmax, lset, title, xvis=False):
    nx, ny = field.shape
    cx, cy = nx/2, ny/2
    dx = minmax[0] - cx
    dy = minmax[1] - cy
    lsetx, lsety = lset.xs, lset.ys
    xs, ys = zip(*region)
    f_cpy = np.zeros_like(field)
    f_cpy[xs, ys] = field[xs, ys]
    f_cpy[lsetx, lsety] = 1.1 * f_cpy.min()
    f_cpy = np.roll(f_cpy, shift=-dx, axis=0)
    f_cpy = np.roll(f_cpy, shift=-dy, axis=1)
    xs, ys = np.where(f_cpy)
    xmin, xmax = xs.min(), xs.max()
    ymin, ymax = ys.min(), ys.max()
    ax.imshow(f_cpy, interpolation='nearest', cmap='hot')
    ax.set_xlim(xmin-2, xmax+2)
    ax.set_ylim(ymin-2, ymax+2)
    ax.set_ylabel(title, size='xx-large', rotation='horizontal')
    if not xvis:
        ax.tick_params(labelbottom=False)

def mag_shear_theta_sectors(psi_arr, cur, den):
    import matplotlib.gridspec as gridspec
    from kaw_analysis.curvature import hessian
    hess = hessian(psi_arr)
    ntheta = 2
    nx, ny = psi_arr.shape
    psi_grad_x, psi_grad_y = gradient(psi_arr)
    bmag = np.sqrt(psi_grad_x**2 + psi_grad_y**2)
    den_x, den_y = gradient(den)
    den_grad_mag = np.sqrt(den_x**2 + den_y**2)
    pl.ion()
    pl.figure()
    pl.imshow(bmag, interpolation='nearest', cmap='hot')
    surf = _cp.TopoSurface(psi_arr)
    regions = surf.get_minmax_regions()
    fig = pl.figure(figsize=(9,12))
    nr = 6
    nc = 2
    for (minmax, pss, type), region in regions.items():
        if len(region) < 100:
            continue
        gs = gridspec.GridSpec(nr, nc, width_ratios=[1,2], hspace=0.0)
        nbr_func = lambda t: _cp.neighbors6(t[0], t[1], nx, ny)
        lset = field_trace._level_set(psi_arr, level_val=psi_arr[pss],
                position=pss, neighbors_func=nbr_func)
        bx = -psi_grad_y
        by = psi_grad_x
        br, btheta = cartesian_to_polar(bx, by, minmax[0], minmax[1])
        region = expand_region_circ(psi_arr, region, minmax, extra=2.0)
        lset = lset.intersection(region)
        # psi
        ax = pl.subplot(gs[0])
        ax.set_title("Field & separatrix", size='x-large')
        field_region_image(ax, psi_arr, region, minmax, lset, title=r'$\psi$')
        # ax = pl.subplot2grid((nr, nc), (0, 1))
        ax = pl.subplot(gs[1])
        ax.set_title(r"Field vs. $r$", size='x-large')
        plot_theta_sectors(ax, psi_arr, minmax, region, ntheta, title=r'$\psi$ vs. $r$')
        # b_theta
        ax = pl.subplot(gs[2])
        field_region_image(ax, btheta, region, minmax, lset, title=r'$B_{\theta}$')
        ax = pl.subplot(gs[3])
        btheta_sectors = plot_theta_sectors(ax, btheta, minmax, region, ntheta, title=r'$B_{\theta}$ vs. $r$')
        # b_shear
        # pl.subplot(nr, nc, 5)
        # field_region_image(btheta, region, minmax, lset, title=r'$B_{\theta}$')
        # pl.subplot(nr, nc, 6)
        # plot_sect_derivs_by_r(btheta_sectors, title=r'$\partial_{r} \frac{B_{\theta}}{r}$ vs. $r$')
        # den
        ax = pl.subplot(gs[4])
        field_region_image(ax, den, region, minmax, lset, title=r'$n$')
        ax = pl.subplot(gs[5])
        plot_theta_sectors(ax, den, minmax, region, ntheta, title=r'$n$ vs. $r$')
        # den_grad
        ax = pl.subplot(gs[6])
        field_region_image(ax, den_grad_mag, region, minmax, lset, title=r'$|\nabla n|$')
        ax = pl.subplot(gs[7])
        plot_theta_sectors(ax, den_grad_mag, minmax, region, ntheta, title=r'$|\nabla n|$ vs. $r$')
        # cur
        ax = pl.subplot(gs[8])
        field_region_image(ax, cur, region, minmax, lset, title=r'$J$')
        ax = pl.subplot(gs[9])
        plot_theta_sectors(ax, cur, minmax, region, ntheta, title=r'$J$ vs. $r$')
        # hessian
        ax = pl.subplot(gs[10])
        field_region_image(ax, hess, region, minmax, lset, title=r'$H(\psi)$', xvis=True)
        ax = pl.subplot(gs[11])
        ax.set_xlabel(r'$r$ (norm. units)', size='x-large')
        plot_theta_sectors(ax, hess, minmax, region, ntheta, title=r'$H(\psi)$ vs. $r$', xvis=True)
        raw_input('enter to continue')
        pl.clf()
    pl.close('all')

def mag_shear_rad_scatter(psi_arr):
    nx, ny = psi_arr.shape
    # mshear = mag_shear(psi_arr)
    psi_grad_x, psi_grad_y = gradient(psi_arr)
    bmag = np.sqrt(psi_grad_x**2 + psi_grad_y**2)
    pl.ion()
    pl.figure()
    pl.imshow(bmag, interpolation='nearest', cmap='hot')
    surf = _cp.TopoSurface(psi_arr)
    regions = surf.get_minmax_regions()
    pl.figure()
    for (minmax, pss, type), region in regions.items():
        br, btheta = cartesian_to_polar(-psi_grad_y, psi_grad_x, minmax[0], minmax[1])
        if len(region) < 100:
            continue
        # region = expand_region(psi_arr, region, ntimes=2)
        region = expand_region_circ(psi_arr, region, minmax, extra=10.0)
        # dists, shears = \
                # flux_tube_radial_scatter(region, minmax, mshear)
        dists, bthetas = \
                flux_tube_radial_scatter(region, minmax, btheta)
        dists, bmags = \
                flux_tube_radial_scatter(region, minmax, bmag)
        # dists, fluxes = \
                # flux_tube_radial_scatter(region, minmax, psi_arr)
        pl.subplot(221)
        pl.scatter(dists, bthetas, c='b', marker='s')
        pl.grid()
        pl.title(r'$B_{\theta}$ vs. $r$')
        pl.subplot(222)
        pl.scatter(dists, bthetas/dists, c='g', marker='d')
        pl.grid()
        pl.title(r'$B_{\theta}/r$ vs. $r$')
        # nonlin_shear = (shears - bmags / dists) / dists
        # pl.scatter(dists, nonlin_shear, c='b', marker='s', label=r'$\nabla_{\perp}B$')
        pl.subplot(224)
        # psi_cpy = np.zeros_like(psi_arr)
        btheta_cpy = np.zeros_like(btheta)
        xs, ys = zip(*region)
        btheta_cpy[xs, ys] = btheta[xs, ys]
        pl.imshow(btheta_cpy, interpolation='nearest', cmap='hot')
        pl.title(r'$B_{\theta}$')
        pl.subplot(223)
        bmag_cpy = np.zeros_like(bmag)
        bmag_cpy[xs, ys] = bmag[xs, ys]
        pl.imshow(bmag_cpy, interpolation='nearest', cmap='hot')
        pl.title(r'$|B|$')
        raw_input('enter to continue')
        pl.clf()
    pl.close('all')

def mag_shear_rad_profs(psi_arr):
    psi_x, psi_y = gradient(psi_arr)
    bmag = np.sqrt(psi_x**2 + psi_y**2)
    psi_arr = gaussian_filter(psi_arr, sigma=2.0, mode='wrap')
    mshear = mag_shear(psi_arr)
    surf = _cp.TopoSurface(psi_arr)
    regions = surf.get_minmax_regions()
    for (minmax, pss, type), region in regions:
        pass
    rprofs = radial_profiles(surf, threshold=25, other_arr=bmag)
    pl.ion()
    pl.figure()
    for minmax, (rprof, region) in rprofs.items():
        # minmax_flux = arr_div[minmax]
        pts, extremum_val, avg_bmags, avg_bmags_errs, avg_dists, avg_dists_errs = \
                zip(*rprof)
        # fluxes = np.abs(np.array(fluxes) - minmax_flux)
        # avg_fluxes = np.abs(np.array(avg_fluxes) - minmax_flux)
        pl.subplot(211)
        pl.plot(avg_dists, avg_bmags, 'd-')
        pl.subplot(212)
        pl.imshow(bmag)
        # pl.plot(avg_dists, avg_shear, 'd-')
        raw_input('enter to continue')
        pl.clf()
    import pdb; pdb.set_trace()
    pl.close('all')

def save_mshear_rad_profs(h5fname):
    dta = tables.openFile(h5fname, 'r')
    for idx, psi_arr in enumerate(dta.walkNodes('/psi', 'Array')):
        if idx < 90: continue
        print psi_arr.name
        psi_arr = psi_arr.read()
        mag_shear_rad_profs(psi_arr)
    dta.close()

if __name__ == '__main__':
    h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10/data.h5'
    # save_mshear_rad_profs(h5fname)
    rad_scatter_driver(h5fname)
    # save_mshear_rad_profs(h5fname)
