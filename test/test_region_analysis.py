import _critical_points as _cp
from region_analysis import radial_profiles

from nose.tools import ok_, eq_

from test_critical_point_network import random_periodic_upsample

def test_radial_profiles():
    arr = random_periodic_upsample(256, 4, seed=0)
    surf = _cp.TopoSurface(arr)
    rprofs = radial_profiles(surf, arr, threshold=50)
    import pdb; pdb.set_trace()
