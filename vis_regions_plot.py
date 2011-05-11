import pylab as pl
from field_trace.region_analysis import expand_regions
from composite_vis import visualize
import _critical_points as _cp
import numpy as np
from itertools import izip
import tables

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
        surf = _cp.TopoSurface(psi_arr)
        regions = surf.get_minmax_regions().values()
        mask = np.zeros_like(bmag)
        mask.fill(-1)
        for region in regions:
            xs, ys = zip(*region)
            mask[xs, ys] = bmag[xs, ys]
        mask[mask==-1] = 1.1 * mask.max()
        pl.figure()
        pl.imshow(bmag, interpolation='nearest', cmap='hot')
        pl.title(r'$|B|$')
        pl.colorbar(pad=0.0, shrink=0.9)
        pl.savefig('%s/bmag-fig.pdf' % savedir)
        visualize(psi_arr, crit_pts=surf.crit_pts,
                 cmap='hot', save_fig=('%s/psi-crit-pts' % savedir),
                 fig_title=r'$\psi$ with critical points', exts=('.pdf',))
        visualize(mask, cmap='hot', save_fig=('%s/bmag-regions-crit-pts' % savedir), exts=('.pdf',),
                 fig_title=r'$|B|$ regions')
    dta.close()

# if __name__ == '__main__':
    # h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10/data.h5'
    # h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10-long-time/data.h5'
    # pl.ioff()
    # vis_regions(h5fname, 2000)
    # save_mshear_rad_profs(h5fname)
