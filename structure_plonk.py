import pylab as pl
import numpy as np
from mag_shear import cartesian_to_polar_coords
from scipy.optimize import leastsq
from scipy.stats import kurtosis

def structure_plonk(size, nstruts, amp, rad_exp, core_radius=5, seed=None):
    arr = np.zeros((size, size), dtype=np.double)
    if seed is not None:
        np.random.seed(seed)
    for _ in range(nstruts):
        # the structure center
        cx, cy = np.random.randint(0, size, 2)
        X, Y, R, theta = cartesian_to_polar_coords(cx, cy, size, size)
        # the distance from the center
        dist = np.sqrt(X**2 + Y**2)
        dist[dist<core_radius] = 0.0
        struct = dist**(-rad_exp)
        struct[np.isinf(struct)] = 0.0
        sign = 1 if np.random.randint(0, 2) else -1
        struct *= (sign * amp / np.max(struct))
        struct[dist<core_radius] = sign * amp
        arr += struct
    arr -= arr.mean()
    return arr

def hists(arrs, labels, title, savename):
    nbins = np.sqrt(arrs[0].size)
    pl.figure()
    fits = [gaus_fit(np.histogram(arr, normed=True, bins=nbins), arr.mean(), arr.std()) \
            for arr in arrs]
    (ns, bins, patches) = pl.hist(x=arrs, histtype='step', label=labels, log=True, bins=100, normed=True, linewidth=3)
    ymin, ymax = pl.ylim()
    pl.ylim(1e-4, ymax)
    for fit in fits:
        best_x ,best_sig = fit
        pl.plot(bins, gaussian(bins, best_x, best_sig), linewidth=3)
    pl.legend()
    pl.grid()
    pl.xlabel('Amplitude (norm. units)')
    pl.ylabel('PDF')
    pl.title(title)
    pl.savefig(savename)

def all_hists(exponent, savename, title, seed=1, repeat=4, size=512):
    np.random.seed(seed)
    exponents = (6, 9)
    arrs = []
    for exp in exponents:
        concat_arrs = []
        for _ in range(repeat):
            concat_arrs.append(structure_plonk(size=size, nstruts=2**exp, amp=1.0, rad_exp=exponent, core_radius=2))
        arrs.append(np.concatenate(concat_arrs).flatten())
    labels = [r'{0:d} strctures, $\kappa={1:3.1f}$'.format(2**i, kurtosis(arr)) for i, arr in zip(exponents, arrs)]
    hists(arrs, labels, title, savename=savename)

def gaussian(Xs, x0, sig):
    norm = 1. / (sig * np.sqrt(2 * np.pi))
    return norm * np.exp(-(Xs-x0)**2 / (2. * sig**2))

def gaus_fit(dta_hist, x0, sig0):
    def resid_func(y):
        Y, X = dta_hist
        x0, sig = y
        expected = gaussian(X[:-1], x0, sig)
        return Y - expected
    y0 = np.array([x0, sig0])
    result, success = leastsq(resid_func, y0)
    if success not in range(1,5):
        raise RuntimeError("result not found")
    return result

if __name__ == '__main__':
    all_hists(exponent=1, repeat=8, savename='rad-exp-1-sim.pdf', title=r'Simulated $r^{-1}$ density structure PDFs')
    all_hists(exponent=2, repeat=8, savename='rad-exp-2-sim.pdf', title=r'Simulated $r^{-2}$ density structure PDFs')
