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

    rt_y = - (a1 + b1 * rt_x) / (c1 + d1 * rt_x)

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

    rt_x = - (a1 + c1 * rt_y) / (b1 + d1 * rt_y)

    if not (0.0 <= rt_x <= 1.0):
        rt_x = - (a2 + c2 * rt_y) / (b2 + d2 * rt_y)

        if not (0.0 <= rt_x <= 1.0):
            raise ValueError("no valid y root found")

    return (rt_x + i, rt_y + j)

def find_cell_zero(deriv, i, j):

    a1, b1, c1, d1 = _get_abcd(deriv.perp_deriv1, i, j)
    a2, b2, c2, d2 = _get_abcd(deriv.perp_deriv2, i, j)

    return _ver1(a1, b1, c1, d1, a2, b2, c2, d2, i, j)

def is_null(xs, ys):
    if _field_trace.same_sign(*xs) or _field_trace.same_sign(*ys):
        return False
    return True

class NullCell(object):

    def __init__(self, deriv, loc):
        self.deriv = deriv
        self._levelset = None
        self._regions = None
        self.loc = loc
        self.bounds = self.deriv.perp_deriv1.shape
        self.x0, self.y0 = find_cell_zero(self.deriv, *self.loc)

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
        if len(self.levelset) <= 3:
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
                    find_cell_zero(deriv, i, j)
                except ValueError:
                    pass
                else:
                    null_cells.append(NullCell(deriv, (i,j)))

    return null_cells

def level_sets(arr, interp, nulls):
    return dict((null, null.levelset) for null in nulls)

def marked_to_mask(shape, marked):
    mask = np.zeros(shape, dtype=np.bool_)
    for m in marked:
        for (i,j) in m:
            mask[i,j] = True
    return mask

def _level_set(arr, level_val, position):
    from collections import deque
    nx, ny = arr.shape
    larr = arr - level_val # leveled array
    i0, j0 = int(position[0]) % nx, int(position[1]) % ny
    cvs = _get_corner_vals(larr, i0, j0)
    # FIXME: this should be the case, is it an interpolation issue?
    # if _field_trace.same_sign_or_zero(*cvs):
        # import pdb; pdb.set_trace()
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

    return np.array(list(marked), dtype=np.int_)

# def neighbors4(i, j, nx, ny):
    # return (((i+1)%nx, j),
            # (i, (j+1)%ny),
            # ((i-1)%nx, j),
            # (i, (j-1)%ny))

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
        regions.append(np.transpose(np.vstack(wh)))
    return regions

def set_region(arr, posns, value):
    posns = np.asarray(posns)
    px = posns[:,0]
    py = posns[:,1]
    arr[px, py] = value

def _is_peak_null(nx, ny, boundary):
    if len(boundary) <= 3:
        return True
    return len(partition_regions(nx, ny, boundary)) == 1

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
