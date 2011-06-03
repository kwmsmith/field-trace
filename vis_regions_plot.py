import pylab as pl
from composite_vis import visualize
import _critical_points as _cp
import numpy as np
from itertools import izip
import tables
from kaw_analysis.curvature import hessian

def mask_arr(arr, regions):
    mask = np.zeros_like(arr)
    mask.fill(-1)
    for region in regions:
        xs, ys = zip(*region)
        mask[xs, ys] = arr[xs, ys]
    mask[mask==-1] = 1.1 * mask.max()
    return mask

def vis_regions(h5fname, idx, savedir='.'):
    idx_str = str(idx)
    dta = tables.openFile(h5fname, 'r')
    arrs = izip(dta.walkNodes('/psi', 'Array'),
                dta.walkNodes('/bx', 'Array'),
                dta.walkNodes('/by', 'Array'))
    for psi_arr, bx_arr, by_arr in arrs:
        if not psi_arr.name.endswith(idx_str):
            continue
        print psi_arr.name
        psi_arr = psi_arr.read()
        psi_arr = psi_arr.astype(np.double)
        bx_arr = bx_arr.read().astype(np.double)
        by_arr = by_arr.read().astype(np.double)
        bmag = bx_arr**2 + by_arr**2
        hess = hessian(psi_arr)
        surf = _cp.TopoSurface(psi_arr)
        regions = surf.get_minmax_regions().values()
        bmag_mask = mask_arr(bmag, regions)
        hess_mask = mask_arr(hess, regions)
        pl.figure()
        pl.imshow(bmag, interpolation='nearest', cmap='hot')
        pl.title(r'$|B|$')
        pl.colorbar(pad=0.0, shrink=0.9)
        pl.savefig('%s/bmag-fig.pdf' % savedir)
        visualize(psi_arr, crit_pts=surf.crit_pts,
                 cmap='hot', save_fig=('%s/psi-crit-pts' % savedir),
                 fig_title=r'$\psi$ with critical points', exts=('.pdf',))
        visualize(bmag_mask, cmap='hot', save_fig=('%s/bmag-regions-crit-pts' % savedir), exts=('.pdf',),
                 fig_title=r'$|B|$ regions')
        visualize(hess, cmap='hot', save_fig=('%s/hessian' % savedir),
                  fig_title=r'$H(\psi)$', exts=('.pdf',))
        hess_neg_mask = hess.copy()
        hess_neg_mask[hess <= 0] = 1.1 * np.min(hess)
        visualize(hess_mask, cmap='hot', save_fig=('%s/hess-regions' % savedir),
                  fig_title=r'$H(\psi)$ regions', exts=('.pdf',))
        visualize(hess_neg_mask, cmap='hot', save_fig=('%s/hess-neg-mask' % savedir),
                  fig_title=r'$H(\psi)$, negative regions masked', exts=('.pdf',))
    dta.close()

# if __name__ == '__main__':
    # h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10/data.h5'
    # h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10-long-time/data.h5'
    # pl.ioff()
    # vis_regions(h5fname, 2000)
    # save_mshear_rad_profs(h5fname)
