from anneal_discrete import Annealer
from contour_tree import wraparound_dist
import numpy as np

def get_dist_map(regions_0, regions_1, nx, ny):
    wdist = wraparound_dist(nx, ny)
    dist_map = np.empty((len(regions_0), len(regions_1)), dtype=np.float_)
    for i, region0 in enumerate(regions_0):
        for j, region1 in enumerate(regions_1):
            dist_map[i, j] = wdist(region0, region1)
    return dist_map

dist_weight = 10.0
area_weight = 1.0

def fitness(reg0, reg1, dist):
    area_diff = abs(reg1[2] - reg0[2])
    return dist_weight * dist + area_weight * area_diff

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

def track_forward_backward(regions_0, regions_1, nx, ny):
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


def track_regions(regions_0, regions_1, nx, ny, nnbrs=5):

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
