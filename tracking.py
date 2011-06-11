import _critical_points as _cp
from scipy.ndimage import gaussian_filter
import tables
from contour_tree import wraparound_dist_vec, wraparound_dist
import numpy as np

from itertools import izip

def get_dist_map(regions_0, regions_1, nx, ny):
    wdist = wraparound_dist_vec(nx, ny)
    dist_map = np.empty((len(regions_0), len(regions_1)), dtype=np.float_)
    x0, y0 = regions_0[:,0], regions_0[:,1]
    x1, y1 = regions_1[:,0], regions_1[:,1]
    for i, region0 in enumerate(regions_0):
        dist_map[i, :] = wdist(region0[0], region0[1], x1, y1)
    return dist_map

dist_weight = 1000.0
area_weight = 0.0
val_weight = 0.0

def fitness(reg0, reg1, dist):
    if reg0.ispeak != reg1.ispeak: return np.inf
    r1area, r0area = reg1.area, reg0.area
    area_diff = abs(r1area - r0area)
    val_diff = abs(reg1.val - reg0.val)
    return dist_weight * dist**3 + area_weight * area_diff + val_diff * val_weight

def get_sim_map(_fitness, dist_map, regions_0, regions_1):
    sim_map = list()
    for i, region0 in enumerate(regions_0):
        sim_map.append(list())
        for j, region1 in enumerate(regions_1):
            sim_map[i].append((j, _fitness(region0, region1, dist_map[i,j])))
    return sim_map

def normalize_regions(all_regions, maxdiff):
    maxdiff = float(maxdiff)
    regions = all_regions[-1]
    max_area = max([r.area for r in regions])
    min_area = min([r.area for r in regions])
    areadiff = np.abs(max_area - min_area)
    area_fac = maxdiff / areadiff
    max_val = max([r.val for r in regions])
    min_val = min([r.val for r in regions])
    valdiff = np.abs(max_val - min_val)
    val_fac = maxdiff / valdiff
    for regions in all_regions:
        for r in regions:
            r.area *= area_fac
            r.val *= val_fac

_extract_regions_cache = {}
def extract_regions(h5fname, sigma=None):
    regions = _extract_regions_cache.get((h5fname, sigma), None)
    if regions is not None:
        return regions
    dta = tables.openFile(h5fname, 'r')
    regions = []
    idx = 0
    for psi_arr, den_arr, cur_arr in izip(
                dta.walkNodes('/psi', 'Array'),
                dta.walkNodes('/den', 'Array'),
                dta.walkNodes('/cur', 'Array')):
        psi_arr_name = psi_arr.name
        psi_arr = psi_arr.read()
        print psi_arr_name
        if sigma is not None:
            psi_arr = gaussian_filter(psi_arr, sigma=sigma, mode='wrap')
        surf = _cp.TopoSurface(psi_arr)
        # surf.simplify_surf_network(surf.get_peak_pit_region_area, threshold=100)
        psi_regions = surf.get_minmax_regions()
        tslice = []
        for (cpt, pss, type), region in psi_regions.items():
            area = len(region)
            minmax = psi_arr[cpt]
            treg = TrackRegion(loc=cpt, area=area, val=minmax, ispeak=(type=='peak'), region=region)
            treg.psi_val = minmax
            treg.den_val = den_arr[cpt]
            treg.cur_val = cur_arr[cpt]
            tslice.append(treg)
        regions.append((idx, tslice))
        idx += 1
    dta.close()
    _extract_regions_cache[h5fname, sigma] = regions
    return regions


class TrackRegion(object):

    def __init__(self, loc, area, val, ispeak, region):
        self.loc = tuple(loc)
        self.area = area
        self.val = val
        self.ispeak = ispeak
        self.region = region

    def __eq__(self, other):
        return (self.loc == other.loc and
                self.area == other.area and
                self.val == other.val and
                self.ispeak == other.ispeak)

    def __hash__(self):
        return hash(self.loc) + hash(self.area) + hash(self.val) + hash(self.ispeak)

def row_to_col_mins(arr, rows, cols):
    assert arr.shape == (len(rows), len(cols))
    dd = {}
    ax_mins = np.min(arr, axis=1)
    minlocs = [np.where(ax_mins[i]==arr[i])[0][0] for i in range(len(ax_mins))]
    for row, ml in zip(rows, minlocs):
        dd[row] = cols[ml]
    return dd

def track_regions_greedy(tslice_to_regions, nx, ny):
    tracks = []
    wdist = wraparound_dist(nx, ny)
    tracks = {}
    # initialize tracks with the 0th tslice
    t0, r0s = tslice_to_regions[0]
    for r0 in r0s:
        tracks[r0] = [(t0, r0)]
    for (t0,r0s), (t1,r1s) in zip(tslice_to_regions, tslice_to_regions[1:]):
        fitness_r0_r1 = []
        for r0 in r0s:
            row = [fitness(r0, r1, wdist(r0.loc, r1.loc)) for r1 in r1s]
            fitness_r0_r1.append(row)
        fitness_r0_r1 = np.array(fitness_r0_r1, dtype=np.double)
        r0_to_r1 = row_to_col_mins(fitness_r0_r1, r0s, r1s)
        r1_to_r0 = row_to_col_mins(fitness_r0_r1.T, r1s, r0s)
        for r0, r1 in r0_to_r1.items():
            if r0 not in tracks:
                import pdb; pdb.set_trace()
            r0_back = r1_to_r0[r1]
            if r0_back is r0:
                # forward <--> back mapping is consistent, add to tracks.
                try:
                    tr = tracks.pop(r0)
                except KeyError:
                    import pdb; pdb.set_trace()
                tr.append((t1, r1))
                tracks[r1] = tr
        for r1 in r1s:
            if r1 not in tracks:
                tracks[r1] = [(t1, r1)]
    return tracks.values()

def chop_tracks(tracks, area_frac=0.1):
    new_tracks = []
    for tr in tracks:
        new_tr = [tr[0]]
        for (idx0, reg0), (idx1, reg1) in zip(tr, tr[1:]):
            if abs(reg0.area-reg1.area) / min(reg0.area, reg1.area) < area_frac:
                new_tr.append((idx1, reg1))
            else:
                new_tracks.append(new_tr)
                new_tr = [(idx1, reg1)]
        new_tracks.append(new_tr)
    return new_tracks

import pylab as pl

def plot_areas_vs_longevity(tracks, nx, ny):
    # wdist = wraparound_dist(nx, ny)
    area_time = []
    for tr in tracks:
        idxs, regs = zip(*tr)
        # mean_area = np.mean([reg.area for reg in regs])
        min_area = min([reg.area for reg in regs])
        # max_area = max([reg.area for reg in regs])
        # gmean_area = gmean([reg.area for reg in regs])
        area_time.append((min_area, len(idxs)))
    areas, times = zip(*area_time)
    pl.scatter(times, areas, c='k', marker='o')
    pl.grid()
    pl.savefig('areas_v_longevity.pdf')

def plot_track_field_vals(tracks, len_cutoff, name_base=''):
    psi_fig, den_fig, cur_fig = pl.figure(), pl.figure(), pl.figure()
    for tr in tracks:
        if len(tr) < len_cutoff: continue
        idxs, regs = zip(*tr)
        pl.figure(psi_fig.number)
        pl.plot(idxs, [reg.psi_val for reg in regs], 's-', hold=True)
        pl.figure(den_fig.number)
        pl.plot(idxs, [reg.den_val for reg in regs], 's-', hold=True)
        pl.figure(cur_fig.number)
        pl.plot(idxs, [reg.cur_val for reg in regs], 's-', hold=True)

    def saver(fig, savename, fieldname):
        pl.figure(fig.number)
        pl.grid()
        xlabel = r'$t$ (norm. units)'
        ylabel = r'{0} at structure core (norm. units)'.format(fieldname)
        title = r'{0} at structure core vs. $t$'.format(fieldname)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.title(title)
        pl.savefig(savename)

    saver(psi_fig, "{0}_psi_v_time.pdf".format(name_base), r'$\psi$')
    saver(den_fig, "{0}_den_v_time.pdf".format(name_base), r'$n$')
    saver(cur_fig, "{0}_cur_v_time.pdf".format(name_base), r'$J$')

def main(h5fname, sigma, savedir='.'):
    nx, ny = 512, 512
    tslice_to_regions = extract_regions(h5fname, sigma=sigma)

    tracks = track_regions_greedy(tslice_to_regions, nx, ny)
    tracks = chop_tracks(tracks, area_frac=0.1)
    plot_areas_vs_longevity(tracks, nx, ny)
    plot_track_field_vals(tracks, len_cutoff=80, name_base='{0}/80'.format(savedir))
    plot_track_field_vals(tracks, len_cutoff=60, name_base='{0}/60'.format(savedir))
    plot_track_field_vals(tracks, len_cutoff=40, name_base='{0}/40'.format(savedir))
    plot_track_field_vals(tracks, len_cutoff=20, name_base='{0}/20'.format(savedir))

if __name__ == '__main__':
    h5fname = '/Users/ksmith/Research/thesis/data/large-rhos2-peak-10/data.h5'
    main(h5fname, sigma=8.0)
