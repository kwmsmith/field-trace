import numpy as np
from mag_shear import cartesian_to_polar_coords

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
    return arr
