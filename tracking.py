from contour_tree import wraparound_dist_vec, wraparound_dist
import numpy as np

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

class TrackRegion(object):

    def __init__(self, loc, area, val, ispeak):
        self.loc = tuple(loc)
        self.area = area
        self.val = val
        self.ispeak = ispeak

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
