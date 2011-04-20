from anneal_discrete import Annealer
from contour_tree import wraparound_dist_vec, wraparound_dist
import numpy as np

from scipy.spatial import cKDTree

def get_dist_map(regions_0, regions_1, nx, ny):
    wdist = wraparound_dist_vec(nx, ny)
    dist_map = np.empty((len(regions_0), len(regions_1)), dtype=np.float_)
    x0, y0 = regions_0[:,0], regions_0[:,1]
    x1, y1 = regions_1[:,0], regions_1[:,1]
    for i, region0 in enumerate(regions_0):
        dist_map[i, :] = wdist(region0[0], region0[1], x1, y1)
    return dist_map

dist_weight = 1.0
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

def sort_sim_map(sim_map):
    for sm in sim_map:
        sm.sort(key=lambda x: x[1])

def _track_forward_backward(regions_0, regions_1, nx, ny):
    dist_map_forward = get_dist_map(regions_0, regions_1, nx, ny)
    dist_map_backward = dist_map_forward.T
    sim_map_forward = get_sim_map(fitness, dist_map_forward, regions_0, regions_1)
    sim_map_backward = get_sim_map(fitness, dist_map_backward, regions_1, regions_0)

    sort_sim_map(sim_map_forward)
    sort_sim_map(sim_map_backward)

    region_map = {}

    for iforward, sm in enumerate(sim_map_forward):
        forward_best, forward_fitness = sm[0]
        backward_best, backward_fitness = sim_map_backward[forward_best][0]
        if backward_best == iforward or forward_fitness <= backward_fitness:
            region_map[iforward] = forward_best
        else:
            raise RuntimeError("finish here!!!")
            # region_map[iforward] = 


    import pdb; pdb.set_trace()

class KDTwraparound(object):

    def __init__(self, pts, nx, ny):
        self.nx, self.ny = nx, ny
        self.all_pts = np.vstack([
            pts+[-nx,-ny], pts+[-nx,0], pts+[-nx,+ny],
            pts+[  0,-ny], pts        , pts+[  0,+ny],
            pts+[+nx,-ny], pts+[+nx,0], pts+[+nx,+ny],
            ])
        self.kdt = cKDTree(self.all_pts)

    def query(self, pts, k):
        d, i = self.kdt.query(pts, k=k, p=2)
        return d, (self.all_pts[i] % [self.nx, self.ny])

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

def track_regions_greedy(all_regions, nx, ny):
    wdist = wraparound_dist(nx, ny)
    normalize_regions(all_regions, nx/2)
    r0s = all_regions[0]
    # area_cutoff = np.median(np.array([r.area for r in r0s]))
    area_cutoff = 0.5
    regions0 = all_regions[0]
    regions0 = sorted(regions0, key=lambda x: x.area)
    regions0 = [r for r in regions0 if r.area > area_cutoff]
    tslice_maps = []
    for idx in range(len(all_regions)-1):
        regions1 = all_regions[idx+1]
        regions1 = sorted(regions1, key=lambda x: x.area)
        regions1 = [r for r in regions1 if r.area > area_cutoff]
        r0_to_r1 = {}
        for region0 in regions0:
            best_r1 = min([(fitness(region0, region1, wdist(region0.loc, region1.loc)), region1)
                for region1 in regions1])[1]
            r0_to_r1[region0] = best_r1
        tslice_maps.append(r0_to_r1)
        regions0 = regions1
    # make tracks out of the tslice_maps list.
    tracks = []
    for idx, map in enumerate(tslice_maps):
        # handle those that are in ``tracks``
        for track in tracks:
            if track[-1][1] in map:
                reg = track[-1][1]
                track.append((idx+1, map[reg]))
                del map[reg]
        # start a new track for everything else
        for reg in map:
            tracks.append([(idx, reg), (idx+1, map[reg])])
    return tracks

def track_regions_anneal(regions_0, regions_1, nx, ny, nnbrs=5):

    dist_map = get_dist_map(regions_0, regions_1, nx, ny)

    sim_map = list()
    for i, region0 in enumerate(regions_0):
        sim_map.append(list())
        for j, region1 in enumerate(regions_1):
            sim_map[i].append((j, fitness(region0, region1, dist_map[i,j])))

    sim_cutoff = min(nnbrs, len(regions_1))

    for sm in sim_map:
        sm.sort(key=lambda x: x[1])
        sm[:] = sm[:sim_cutoff]

    tot_area = regions_1[:,2].sum()

    def energy(region_map):
        E = 0.0e0
        seen = set()
        marked_area = 0
        for idx1, idx2 in enumerate(region_map):
            E += fitness(regions_0[idx1], regions_1[idx2], dist_map[idx1, idx2])
            if idx2 in seen:
                continue
            seen.add(idx2)
            marked_area += regions_1[idx2][2]
        assert len(seen) <= len(regions_1)
        assert marked_area <= tot_area
        E += area_weight * (tot_area - marked_area)
        E += len(regions_1) - len(seen)
        return E

    def move(region_map):
        from random import randrange
        for _ in range(1):
            idx = randrange(len(region_map))
            swap_idx = randrange(sim_cutoff)
            region_map[idx] = sim_map[idx][swap_idx][0]

    initial_state = np.array([sm[0][0] for sm in sim_map])

    dEmax = 0.0
    state_cpy = initial_state.copy()
    prev_E = energy(state_cpy)
    for i in range(100):
        move(state_cpy)
        dE = abs(energy(state_cpy) - prev_E)
        if dE > dEmax:
            dEmax = dE

    Tmax = dEmax

    an = Annealer(energy, move)

    state_cpy = initial_state.copy()
    print "initial energy: %f" % energy(initial_state)
    best_state, best_energy = an.anneal_nr(
                                    state=initial_state,
                                    Tmax=Tmax,
                                    temp_steps=50,
                                    Tfac=0.5,
                                    nover=100,
                                    # Tmin=Tmin,
                                    # steps=10**4,
                                    updates=True)
    print "best energy: %f" % best_energy
    print "initial energy: %f" % energy(state_cpy)
    best_state, best_energy = an.anneal(
                                    state=state_cpy,
                                    Tmax=Tmax,
                                    Tmin=1e-12,
                                    steps=5e3,
                                    updates=100)
    print "best energy: %f" % best_energy

    import pdb; pdb.set_trace()
