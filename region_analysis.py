import _critical_points as _cp
import pylab as pl
import numpy as np

from contour_tree import wraparound_dist_vec
from field_trace import Region

def crit_pt_fitness(cpts, mod_grad):
    cpt2grad = {}
    for cpt in cpts:
        cpt2grad[cpt] = mod_grad[cpt]
    return cpt2grad

def region_areas_hist(surf):
    cpts2regions = surf.get_minmax_regions()
    cpts2areas = dict((c, len(r)) for (c,r) in cpts2regions.items())
    nbins = sqrt(len(cpts2areas))
    pl.hist(cpts2areas.values(), bins=nbins, normed=True, histtype='step', linewidth=2)

def region_energy(regionx, regiony, arr):
    return np.sum(arr[regionx, regiony]**2)

def energy_of_regions(regions, arr):
    r2es = []
    for region in regions:
        regionx, regiony = zip(*region)
        r2es.append((region, region_energy(regionx, regiony, arr)))
    return r2es

def region_max(regionx, regiony, arr):
    arr_region = arr[regionx, regiony]
    rmax = arr_region.max()
    rmin = arr_region.min()
    return max(abs(rmax), abs(rmin))

def region_avg(regionx, regiony, arr):
    arr_region = arr[regionx, regiony]
    ave = arr_region.mean()
    return ave

def region_pseudo_volume(regionx, regiony, arr):
    rmax = region_max(regionx, regiony, arr)
    return rmax * len(regionx)

def get_seed_points(region, arr, nseeds, reverse=False):
    if nseeds <= 0:
        nseeds = 1
    ordered_pts = sorted([(arr[pt], pt) for pt in region], reverse=reverse)
    step = max(int(len(ordered_pts)/nseeds), 1)
    # ensure that we include the first & last points...
    # seeds = [ordered_pts[0]] + ordered_pts[1:-2:step] + [ordered_pts[-1]]
    seeds = ordered_pts[0:-1:step]
    return [s[1] for s in seeds]

def flux_tube_radial_scatter(region, minmax_pt, arr):
    wdist = wraparound_dist_vec(*arr.shape)
    region = list(region)
    dists = []
    vals = []
    region_arr = np.array(region, dtype=np.int)
    xs = region_arr[:,0]
    ys = region_arr[:,1]
    dists = wdist(minmax_pt[0], minmax_pt[1], xs, ys)
    vals = arr[xs, ys]
    return (dists, vals)

def flux_tube_radial_spokes(region, minmax_pt, arr, peak_or_pit):
    nx, ny = arr.shape
    ctr_max = nx+ny
    frontier = set()
    for pt in region:
        for nbr in _cp.neighbors6(pt[0], pt[1], nx, ny):
            if nbr not in region:
                frontier.add(pt)
    spokes = []
    for pt in frontier:
        cur = pt
        spoke = [cur]
        if peak_or_pit == 'peak': next = _cp._get_highest_neighbor
        elif peak_or_pit == 'pit': next = _cp._get_lowest_neighbor
        else: raise ValueError("peak_or_pit not in ('peak', 'pit')")
        ctr = ctr_max
        while ctr > 0:
            cur = next(arr, cur)
            spoke.append(cur)
            if cur == minmax_pt:
                spokes.append(spoke)
                break
            ctr -= 1
    return spokes

def flux_tube_radial_profile(region, minmax_pt, surf, peak_or_pit,
        other_arr=None, mask=None):
    arr = surf.arr
    wdist = wraparound_dist_vec(*arr.shape)
    from field_trace import _level_set
    if other_arr is None:
        other_arr = arr
    seed_points = get_seed_points(region, arr, int(np.sqrt(len(region))),
            reverse=(peak_or_pit == 'peak'))
    radial_profile = []
    nx, ny = arr.shape
    for seed in seed_points:
        seed_flux = arr[seed]
        nbr_func = lambda t: _cp.neighbors6(t[0], t[1], nx, ny)
        lset = _level_set(arr, level_val=seed_flux,
                position=seed, neighbors_func=nbr_func)
        lset = lset.intersection(region)
        if not lset:
            continue
        if mask is not None:
            mask[lset.xs, lset.ys] = True
        lset_arr = other_arr[lset.xs, lset.ys]
        lset_mean = lset_arr.mean()
        lset_err = lset_arr.std() / np.sqrt(len(lset_arr))
        lset_dists = wdist(lset.xs, lset.ys, minmax_pt[0], minmax_pt[1])
        dists_mean = lset_dists.mean()
        dists_err = lset_dists.std() / np.sqrt(len(lset_dists))
        radial_profile.append(
                (seed, seed_flux, lset_mean, lset_err, dists_mean, dists_err))
    return radial_profile

def expand_region_circ(arr, region, minmax, extra=0.0):
    wdist = wraparound_dist_vec(*arr.shape)
    region_arr = np.array(list(region), dtype=np.int)
    xs = region_arr[:,0]
    ys = region_arr[:,1]
    dists = wdist(minmax[0], minmax[1], xs, ys)
    maxdist = max(dists)
    return expand_region_to_radius(arr, region, maxdist+extra, minmax)

def expand_region_to_radius(arr, region, radius, minmax):
    nx, ny = arr.shape
    wdist = wraparound_dist_vec(nx, ny)
    frontier = set()
    new_region = set(region)
    for pt in region:
        nbrs = _cp.neighbors6(pt[0], pt[1], nx, ny)
        for nbr in nbrs:
            if nbr not in region and wdist(minmax[0], minmax[1], nbr[0], nbr[1]) <= radius:
                frontier.add(nbr)
                new_region.add(nbr)
    # new_region.update(frontier)
    while frontier:
        new_frontier = set()
        for pt in frontier:
            nbrs = _cp.neighbors6(pt[0], pt[1], nx, ny)
            for nbr in nbrs:
                if nbr not in new_region and wdist(minmax[0], minmax[1], nbr[0], nbr[1]) <= radius:
                    new_frontier.add(nbr)
        frontier = new_frontier
        new_region.update(frontier)
    return new_region

def expand_region(arr, region, ntimes):
    nx, ny = arr.shape
    frontier = set()
    new_region = set(region)
    for pt in region:
        nbrs = _cp.neighbors6(pt[0], pt[1], nx, ny)
        for nbr in nbrs:
            if nbr not in region:
                frontier.add(nbr)
    new_region.update(frontier)
    for _ in range(ntimes-1):
        new_frontier = set()
        for pt in frontier:
            nbrs = _cp.neighbors6(pt[0], pt[1], nx, ny)
            for nbr in nbrs:
                if nbr not in new_region:
                    new_frontier.add(nbr)
        frontier = new_frontier
        new_region.update(frontier)
    return new_region

def radial_profiles(surf, threshold, other_arr=None, expand_regions=0, mask=None):
    regions = surf.get_minmax_regions()
    cpt2rprof = {}
    for ((cpt, peak_or_pit), region) in regions.items():
        if len(region) < threshold:
            continue
        if expand_regions:
            region = expand_region(surf.arr, region, expand_regions)
        rprofile = flux_tube_radial_profile(region, cpt, surf,
                peak_or_pit=peak_or_pit, other_arr=other_arr, mask=mask)
        cpt2rprof[cpt] = (rprofile, region)
    return cpt2rprof
