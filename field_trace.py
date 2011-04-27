import time
from collections import defaultdict

import numpy as np
from numpy.linalg import eig
from numpy.linalg import LinAlgError

from kaw_analysis import vcalc
from wrap_gsl_interp2d import Interp2DPeriodic


import _field_trace

class Derivator(object):

    def __init__(self, scalar_arr):
        self.arr = scalar_arr
        len1 = self.len1 = float(self.arr.shape[0])
        len2 = self.len2 = float(self.arr.shape[1])

        self.arr_interp = Interp2DPeriodic(len1, len2, self.arr)

        self.deriv1 = vcalc.cderivative(self.arr, 'X_DIR')
        self.deriv2 = vcalc.cderivative(self.arr, 'Y_DIR')

        self.perp_deriv1 =  self.deriv2
        self.perp_deriv2 = -self.deriv1

        self.perp_deriv1_interp = Interp2DPeriodic(len1, len2, self.perp_deriv1)
        self.perp_deriv2_interp = Interp2DPeriodic(len1, len2, self.perp_deriv2)

        self.deriv12 = vcalc.cderivative(self.deriv1, 'Y_DIR')
        self.deriv11 = vcalc.cderivative(self.arr, 'X_DIR', order=2)
        self.deriv22 = vcalc.cderivative(self.arr, 'Y_DIR', order=2)

        self.deriv12_interp = Interp2DPeriodic(len1, len2, self.deriv12)
        self.deriv11_interp = Interp2DPeriodic(len1, len2, self.deriv11)
        self.deriv22_interp = Interp2DPeriodic(len1, len2, self.deriv22)

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

def grad_zero_broyden(dd, xin, yin):
    from scipy.optimize import broyden3
    def __eval__(X):
        return (dd.perp_deriv1_interp.eval(*X), dd.perp_deriv2_interp.eval(*X))
    return broyden3(__eval__, (xin, yin), iter=10)

def hessian_esys(derivs, x0, y0):
    val_12 = derivs.deriv12_interp.eval(x0, y0)
    mm = [[derivs.deriv11_interp.eval(x0, y0), val_12],
          [val_12, derivs.deriv22_interp.eval(x0, y0)]]
    evals, evecs = eig(mm)
    return evals, evecs
                   

class NullCell(object):

    def __init__(self, deriv, loc, locmin=None):
        self.deriv = deriv
        self.loc = loc
        self.bounds = self.deriv.perp_deriv1.shape
        if locmin is None:
            self.x0, self.y0 = grad_zero_broyden(self.deriv, *self.loc)
        else:
            self.x0, self.y0 = locmin
        self.esys = hessian_esys(self.deriv, self.x0, self.y0)
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

    def is_maximum(self):
        evs, evecs = self.esys
        return evs[0] < 0 and evs[1] < 0

    def is_minimum(self):
        evs, evecs = self.esys
        return evs[0] > 0 and evs[1] > 0

    def is_saddle(self):
        evs, evecs = self.esys
        return evs[0] * evs[1] < 0

_get_corner_vals = _field_trace._get_corner_vals

def _get_abcd(f, i, j):
    f00, f01, f10, f11 = _get_corner_vals(f, i, j)

    a = f00
    b = f10 - f00
    c = f01 - f00
    d = f11 - f10 - f01 + f00

    return (a, b, c, d)

def find_null_cells_minimize(deriv, thresh_frac=0.05):
    from scipy.optimize import fmin_powell
    grad2 = deriv.deriv1**2 + deriv.deriv2**2
    grad2_interp = Interp2DPeriodic(grad2.shape[0], grad2.shape[1], grad2)
    nx, ny = grad2.shape
    cx, cy = np.where(grad2 < thresh_frac * grad2.max())
    null_cells = []
    marked_cells = set()
    def __eval__(X):
        return grad2_interp.eval(*X)
    for (i,j) in zip(cx, cy):
        try:
            xmin, ymin = fmin_powell(__eval__, (i,j), disp=0)
        except LinAlgError:
            continue
        if (int(xmin), int(ymin)) not in marked_cells:
            try:
                nc = NullCell(deriv, (i,j), locmin=(xmin, ymin))
            except LinAlgError:
                continue
            null_cells.append(nc)
            marked_cells.add((int(xmin), int(ymin)))
    return null_cells

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
                except LinAlgError:
                    pass
    return null_cells

def marked_to_mask(shape, marked):
    mask = np.zeros(shape, dtype=np.bool_)
    for m in marked:
        mask[m.xs, m.ys] = True
    return mask

def _level_set(arr, level_val, position, neighbors_func=None):
    nx, ny = arr.shape
    if neighbors_func is None:
        neighbors_func = lambda t: _field_trace.neighbors8(t[0], t[1], nx, ny)
    i0, j0 = int(position[0]) % nx, int(position[1]) % ny
    marked = set([(i0, j0)])
    front = [(i0, j0)]
    while front:
        pt = front.pop()
        for (i, j) in neighbors_func(pt):
        # for (i, j) in _field_trace.neighbors8(i0, j0, nx, ny):
            if (i,j) in marked:
                continue
            try:
                cvs = _get_corner_vals(arr, i, j)
            except OverflowError:
                print i, j
                raise
            if not _field_trace.same_sign_or_zero(*cvs, lval=level_val):
                marked.add((i, j))
                front.append((i, j))

    psns = np.array(list(marked))
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
        self.xs = np.asanyarray(xs, dtype=np.int64)
        self.ys = np.asanyarray(ys, dtype=np.int64)
        self.ncontained = 0
        self.contained_nulls = set()
        self._region_set = None

    def _get_region_set(self):
        rset = set()
        for x, y in zip(self.xs, self.ys):
            rset.add((x,y))
        return rset

    region_set = property(_get_region_set)

    def is_subregion(self, other):
        return self.size <= other.size and \
                np.any((self.xs[0] == other.xs) & (self.ys[0] == other.ys))

    def __contains__(self, (x, y)):
        return (x, y) in self.region_set

    def intersection(self, other):
        new_set = self.region_set.intersection(other)
        xs, ys = zip(*new_set)
        return Region(xs, ys)


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

def find_region_ncontained(shape, regions):
    # '''
    # given a collection `regions`, finds the collection of all regions `rset`
    # contained in each region `R` and sets `R.contains` to `rset`.

    # '''
    NULL_FLAG = -1
    map_arr = np.zeros(shape, dtype=np.int32)
    map_arr.fill(NULL_FLAG)
    regions = regions[:]
    regions.sort(key=lambda r: r.size)
    for idx, region in enumerate(regions):
        contained_idxs = np.unique(map_arr[region.xs, region.ys])
        # `contained_idxs` is a *sorted* array of indices into `regions`.
        for cidx in contained_idxs:
            if cidx == NULL_FLAG:
                continue
            assert cidx >= 0
            region.ncontained += regions[cidx].ncontained + 1
            # region.contains.append(regions[cidx])
        map_arr[region.xs, region.ys] = idx

def null_cover_regions(nulls, regions):
    '''
    Filters `regions` to the smallest set of regions that contain all nulls in
    `nulls`.

    '''
    EMPTY_VAL = -1
    shape = nulls[0].bounds
    null_mask = np.empty(shape, dtype=np.int64)
    null_mask.fill(EMPTY_VAL)
    uncovered_nulls = set(range(len(nulls)))
    regions = regions[:]
    regions.sort(key=lambda r: r.size, reverse=True)
    covering_regions = []
    for idx, null in enumerate(nulls):
        null_mask[null.loc[0], null.loc[1]] = idx
    while uncovered_nulls:
        region = regions.pop()
        contained_nulls = np.unique(null_mask[region.xs, region.ys])
        for null_idx in contained_nulls:
            if null_idx in uncovered_nulls:
                uncovered_nulls.remove(null_idx)
                region.contained_nulls.add(null_idx)
        if region.contained_nulls:
            covering_regions.append(region)
    return covering_regions

def filter_min_regions(shape, regions, min_size=0):
    '''
    given a collection of regions, returns a list of all minimal regions,
    defined to be the set of regions that contain no other regions.

    All regions smaller than min_size are ignored.

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

def get_nulls(arr):
    arr = np.asanyarray(arr, dtype=np.double)
    dd = Derivator(arr)
    # return find_null_cells(dd)
    return find_null_cells_minimize(dd, thresh_frac=0.02)

def nulls_and_regions(arr, chatty=False):
    arr = np.asanyarray(arr, dtype=np.double)
    dd = Derivator(arr)
    if chatty:
        print "%s: locating nulls" % (time.asctime())
    nulls = find_null_cells(dd)
    if chatty:
        print "%s: number of nulls: %d" % (time.asctime(), len(nulls))
    if chatty:
        print "%s: getting regions" % (time.asctime(),)
    regions = []
    for null in nulls:
        regions.extend(null.regions)
    if chatty:
        print "%s: number of regions: %d" % (time.asctime(), len(regions))
    return (nulls, regions)

def regions_by_n_contained(regions):
    n_to_regions = defaultdict(set)
    for region in regions:
        n_to_regions[region.ncontained].add(region)
        # n_to_regions[len(region.contains)].add(region)
    return n_to_regions

def detect_min_regions(arr, min_size=0):
    regions = all_regions(arr, chatty=True)

    min_regions = filter_min_regions(arr.shape, regions, min_size=min_size)

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

def save_fig_with_scatter(arr, scatter_xy, basename, fig_exts=('.eps', '.png', '.pdf')):
    import pylab as pl
    pl.ioff()
    fig = pl.figure()
    pl.imshow(arr, interpolation='nearest', cmap='hot')
    X, Y = scatter_xy
    pl.scatter(Y, X, c='b')
    for ext in fig_exts:
        pl.savefig(basename+ext)
    pl.close(fig)
