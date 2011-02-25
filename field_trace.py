from collections import defaultdict

import numpy as np

from kaw_analysis import vcalc
from wrap_gsl_interp2d import Interp2DPeriodic

from nose.tools import set_trace

import _field_trace

class Derivator(object):

    def __init__(self, scalar_arr, len1, len2):
        self.arr = scalar_arr
        self.len1 = float(len1)
        self.len2 = float(len2)

        self.arr_interp = Interp2DPeriodic(len1, len2, self.arr)

        self.deriv1 = vcalc.cderivative(self.arr, 'X_DIR')
        self.deriv2 = vcalc.cderivative(self.arr, 'Y_DIR')

        self.perp_deriv1 =  self.deriv2
        self.perp_deriv2 = -self.deriv1

        self.perp_deriv1_interp = Interp2DPeriodic(len1, len2, self.perp_deriv1)
        self.perp_deriv2_interp = Interp2DPeriodic(len1, len2, self.perp_deriv2)

        deriv12 = vcalc.cderivative(self.deriv1, 'Y_DIR')
        deriv11 = vcalc.cderivative(self.arr, 'X_DIR', order=2)
        deriv22 = vcalc.cderivative(self.arr, 'Y_DIR', order=2)

        self.deriv12_interp = Interp2DPeriodic(len1, len2, deriv12)
        self.deriv11_interp = Interp2DPeriodic(len1, len2, deriv11)
        self.deriv22_interp = Interp2DPeriodic(len1, len2, deriv22)

def quad_roots(A, B, C):
    '''solve quadratic equation for A x**2 + B x + C == 0'''
    from math import sqrt
    B /= A
    C /= A

    discr = sqrt(B**2 - 4.0 * C)
    if B < 0.0:
        rt = -B + discr
    elif B > 0.0:
        rt = -B - discr

    rt /= 2.0

    other_rt = C / rt

    return (rt, other_rt)

def _ver0(a1, b1, c1, d1, a2, b2, c2, d2, i, j):
    C = a1 * c2 - c1 * a2
    B = (a1 * d2 - d1 * a2) + (b1 * c2 - c1 * b2)
    A = b1 * d2 - d1 * b2

    rt_x0, rt_x1 = quad_roots(A, B, C)

    if 0.0 <= rt_x0 <= 1.0:
        rt_x = rt_x0
    elif 0.0 <= rt_x1 <= 1.0:
        rt_x = rt_x1
    else:
        raise ValueError("no valid x root found")

    try:
        rt_y = - (a1 + b1 * rt_x) / (c1 + d1 * rt_x)
    except ZeroDivisionError:
        rt_y = -1.0

    if not (0.0 <= rt_y <= 1.0):
        rt_y = - (a2 + b2 * rt_x) / (c2 + d2 * rt_x)

        if not (0.0 <= rt_y <= 1.0):
            raise ValueError("no valid y root found")

    return (rt_x + i, rt_y + j)

def _ver1(a1, b1, c1, d1, a2, b2, c2, d2, i, j):
    C = a1 * b2 - b1 * a2
    B = (a1 * d2 - d1 * a2) - (b1 * c2 - c1 * b2)
    A = c1 * d2 - d1 * c2

    rt_y0, rt_y1 = quad_roots(A, B, C)

    if 0.0 <= rt_y0 <= 1.0:
        rt_y = rt_y0
    elif 0.0 <= rt_y1 <= 1.0:
        rt_y = rt_y1
    else:
        raise ValueError("no valid y root found")

    try:
        rt_x = - (a1 + c1 * rt_y) / (b1 + d1 * rt_y)
    except ZeroDivisionError:
        rt_x = -1.0

    if not (0.0 <= rt_x <= 1.0):
        rt_x = - (a2 + c2 * rt_y) / (b2 + d2 * rt_y)

        if not (0.0 <= rt_x <= 1.0):
            raise ValueError("no valid y root found")

    return (rt_x + i, rt_y + j)

def find_cell_zero(deriv, i, j):

    a1, b1, c1, d1 = _get_abcd(deriv.perp_deriv1, i, j)
    a2, b2, c2, d2 = _get_abcd(deriv.perp_deriv2, i, j)

    try:
        return _ver1(a1, b1, c1, d1, a2, b2, c2, d2, i, j)
    except ZeroDivisionError:
        return _ver0(a1, b1, c1, d1, a2, b2, c2, d2, i, j)

def is_null(xs, ys):
    if _field_trace.same_sign(*xs) or _field_trace.same_sign(*ys):
        return False
    return True

class NullCell(object):

    def __init__(self, deriv, loc):
        self.deriv = deriv
        self.loc = loc
        self.bounds = self.deriv.perp_deriv1.shape
        self.x0, self.y0 = find_cell_zero(self.deriv, *self.loc)
        self._levelset = None
        self._regions = None

    def get_levelset(self):
        if self._levelset is None:
            level_value = self.deriv.arr_interp.eval(self.x0, self.y0)
            self._levelset = _level_set(self.deriv.arr, level_value, (self.x0, self.y0))
        return self._levelset

    levelset = property(get_levelset)

    def get_regions(self):
        if self._regions is None:
            self._regions = partition_regions(self.bounds[0], self.bounds[1], self.levelset)
        return self._regions

    regions = property(get_regions)

    def is_saddle(self):
        if self.levelset.size <= 3:
            return False
        return len(self.regions) > 1

_get_corner_vals = _field_trace._get_corner_vals

def _get_abcd(f, i, j):
    f00, f01, f10, f11 = _get_corner_vals(f, i, j)

    a = f00
    b = f10 - f00
    c = f01 - f00
    d = f11 - f10 - f01 + f00

    return (a, b, c, d)

def find_null_cells(deriv):
    null_cells = []
    nx, ny = deriv.perp_deriv1.shape
    assert (nx, ny) == deriv.perp_deriv2.shape
    for i in range(nx):
        for j in range(ny):
            if is_null(_get_corner_vals(deriv.perp_deriv1, i, j),
                       _get_corner_vals(deriv.perp_deriv2, i, j)):
                try:
                    null_cells.append(NullCell(deriv, (i,j)))
                except ValueError:
                    pass

    return null_cells

def marked_to_mask(shape, marked):
    mask = np.zeros(shape, dtype=np.bool_)
    for m in marked:
        mask[m.xs, m.ys] = True
    return mask

def _level_set(arr, level_val, position):
    from collections import deque
    nx, ny = arr.shape
    larr = arr - level_val # leveled array
    i0, j0 = int(position[0]) % nx, int(position[1]) % ny
    cvs = _get_corner_vals(larr, i0, j0)
    marked = set([(i0, j0)])
    front = [(i0, j0)]
    while front:
        i0, j0 = front.pop()
        for (i, j) in _field_trace.neighbors8(i0, j0, nx, ny):
            if (i,j) in marked:
                continue
            cvs = _get_corner_vals(larr, i, j)
            if not _field_trace.same_sign_or_zero(*cvs):
                marked.add((i, j))
                front.append((i, j))

    psns = np.array(list(marked), dtype=np.uint16)
    return Region(psns[:,0], psns[:,1])

neighbors4 = _field_trace.neighbors4

def partition_regions(nx, ny, boundary):
    '''
    returns a list of np arrays of indicies.  Each array of indices is a
    connected region broken up by the boundary.

    '''
    EMPTY_COLOR = 1
    BOUNDARY_COLOR = 0
    domain = np.zeros((nx, ny), dtype=np.int32)
    domain.fill(EMPTY_COLOR)
    set_region(domain, boundary, BOUNDARY_COLOR)
    larr = connected_component_label(domain)
    maxlabel = larr.max()
    regions = []
    for label in range(1, maxlabel+1):
        wh = np.where(larr == label)
        if len(wh[0]):
            regions.append(Region(*wh))
    return regions

class Region(object):

    def __init__(self, xs, ys):
        if len(xs) != len(ys):
            raise ValueError("invalid coordinates")
        self.size = len(xs)
        self.xs = np.asanyarray(xs, dtype=np.uint16)
        self.ys = np.asanyarray(ys, dtype=np.uint16)

    def is_subregion(self, other):
        return self.size <= other.size and \
                np.any((self.xs[0] == other.xs) & (self.ys[0] == other.ys))

    '''
    def _is_subregion(self, other):
        return self.size <= other.size and \
                np.any(np.intersect1d(np.where(self.xs[0]==other.xs)[0],
                                      np.where(self.ys[0]==other.ys)[0]))


    def _is_subregion(self, other):
        if self.size > other.size:
            return False
        idxs = np.transpose(np.vstack([other.xs, other.ys]))
        return [self.xs[0], self.ys[0]] in idxs
    '''

def set_region(arr, region, value):
    arr[region.xs, region.ys] = value

def connected_component_label(arr, output=None):
    '''
    arr -- array object; non-zero values are treated as features and zero
    values are treated as background.

    if output is supplied, will set output's values such that connected
    components in arr are labeled with a single value in output.

    Arr is considered to be a periodic array, so edges connect.

    '''
    from scipy.ndimage import label
    from collections import defaultdict

    if output is not None:
        larr = output
        nlabels = label(arr, output=larr)
    else:
        larr, nlabels = label(arr, output=output)

    equiv_classes = dict((v, set([v])) for v in range(1, nlabels+1))

    compares = ((larr[0,:], larr[-1,:]), # top & bottom rows
                (larr[:,0], larr[:,-1])) # first & last columns

    for cs in compares:
        for v0, v1 in zip(*cs):
            if v0 and v1 and v0 != v1:
                union = equiv_classes[v0].union(equiv_classes[v1])
                equiv_classes[v0] = equiv_classes[v1] = union

    marked = set()
    for ec in equiv_classes:
        cl = equiv_classes[ec]
        if len(cl) > 1:
            minval = min(cl)
            for v in cl:
                if v != minval and v not in marked:
                    larr[larr == v] = minval
                    marked.add(v)

    return larr

# def _filter_min_regions(regions, min_size=0):
    # '''
    # given an iterable of regions, returns a list of all minimal regions,
    # defined to be the set of regions that contain no other regions.

    # All regions smaller than min_size are treated as though they don't exist.

    # '''
    # regions = regions[:]
    # min_regions = []
    # regions.sort(key=lambda r: r.size, reverse=True)
    # while regions:
        # cur_min = regions.pop()
        # if cur_min.size < min_size:
            # continue
        # min_regions.append(cur_min)
        # new_regions = []
        # for region in regions:
            # if not cur_min.is_subregion(region):
                # new_regions.append(region)
        # regions = new_regions
    # return min_regions

def filter_min_regions(shape, regions, min_size=0):
    '''
    given an iterable of regions, returns a list of all minimal regions,
    defined to be the set of regions that contain no other regions.

    All regions smaller than min_size are treated as though they don't exist.

    '''
    mask = np.zeros(shape, dtype=np.bool_)
    regions = regions[:]
    min_regions = []
    regions.sort(key=lambda r: r.size)
    for region in regions:
        if region.size < min_size:
            continue
        if not np.any(mask[region.xs, region.ys]):
            mask[region.xs, region.ys] = True
            min_regions.append(region)
    return min_regions

def detect_min_regions(arr, min_size=0):
    import time
    arr = np.asanyarray(arr, dtype=np.double)
    N = arr.shape[0]

    print "%s: locating nulls" % (time.asctime())
    dd = Derivator(arr, N, N)
    nulls = find_null_cells(dd)
    print "%s: number of nulls: %d" % (time.asctime(), len(nulls))

    print "%s: getting min regions" % (time.asctime())
    regions = []
    for null in nulls:
        regions.extend(null.regions)

    print "%s: number of regions: %d" % (time.asctime(), len(regions))

    min_regions = filter_min_regions(arr.shape, regions, min_size=min_size)

    print "%s: number of min regions: %d" % (time.asctime(), len(min_regions))

    return min_regions

def regions_to_mask(shape, regions):
    mask = np.zeros(shape, dtype=np.bool_)
    for region in regions:
        mask[region.xs, region.ys] = True
    return mask

def save_fig(arr, basename, fig_exts=('.eps', '.png', '.pdf')):
    import pylab as pl
    pl.ioff()
    fig = pl.figure()
    pl.imshow(arr, interpolation='nearest', cmap='hot')
    for ext in fig_exts:
        pl.savefig(basename+ext)
    pl.close(fig)
