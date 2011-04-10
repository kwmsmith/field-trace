import pylab as pl
import numpy as np

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

def get_seed_points(region, arr, nseeds):
    ordered_pts = sorted([(arr[pt], pt) for pt in region])
    step = max(int(len(ordered_pts)/nseeds), 1)
    # ensure that we include the first & last points...
    seeds = [ordered_pts[0]] + ordered_pts[1:-2:step] + [ordered_pts[-1]]
    return [s[1] for s in seeds]
    
def flux_tube_radial_profile(region, min_max_pt, arr, other_arr=None):
    from field_trace import _level_set
    other_arr = other_arr or arr
    seed_points = get_seed_points(region, arr, int(np.sqrt(len(region)))+5)
    radial_profile = []
    for seed in seed_points:
        seed_flux = arr[seed]
        lset = _level_set(arr, level_val=seed_flux, position=seed)
        lset_avg = other_arr[lset.xs, lset.ys].mean()
        radial_profile.append((seed, seed_flux, lset_avg))
    return radial_profile

def radial_profiles(surf, arr, threshold, other_arr=None):
    cpts2regions = surf.get_minmax_regions()
    cpt2rprof = {}
    for cpt, region in cpts2regions.items():
        if len(region) < threshold:
            continue
        rprofile = flux_tube_radial_profile(region, cpt, arr, other_arr=other_arr)
        cpt2rprof[cpt] = rprofile
    return cpt2rprof
