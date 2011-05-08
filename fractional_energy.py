import numpy as np
import tables
from tracking import extract_regions

from test_tracking import h5fname

def fractional_energy(h5fname, sigma):
    tslice_to_regions = extract_regions(h5fname, sigma=sigma)
    dta = tables.openFile(h5fname, 'r')
    efrac = []
    rfrac = []
    # efrac_by_regfrac = []
    # dual = []
    bx_arrs = dta.walkNodes('/bx', 'Array')
    by_arrs = dta.walkNodes('/by', 'Array')
    for idx, (bx_arr, by_arr) in enumerate(zip(bx_arrs, by_arrs)):
        bx_arr = bx_arr.read()
        by_arr = by_arr.read()
        bmag_arr = bx_arr**2 + by_arr**2
        bmag_arr_regions = np.zeros(bx_arr.shape, dtype=np.bool)
        for reg in tslice_to_regions[idx][1]:
            bmag_arr_regions[zip(*reg.region)] = True
        reg_frac = float(np.sum(bmag_arr_regions)) / bmag_arr_regions.size
        eng_frac = np.sum(bmag_arr[bmag_arr_regions]) / np.sum(bmag_arr)
        efrac.append(eng_frac)
        rfrac.append(reg_frac)
    dta.close()
    return np.array(efrac), np.array(rfrac)

import pylab as pl
pl.ion()
efrac, rfrac = fractional_energy('./data.h5', sigma=8)
pl.plot(efrac, 'ro-', label=r'$\Delta E$')
pl.plot(rfrac, 'kd-', label=r'$\Delta A$', hold=True)
pl.plot(efrac/rfrac, 'gs-', label=r'$\frac{\Delta E}{\Delta A}$', hold=True)
pl.grid()
pl.legend()
raw_input('enter to continue')
