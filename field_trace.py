from collections import defaultdict

import numpy as np

from kaw_analysis import vcalc
from wrap_gsl_interp2d import Interp2DPeriodic

from nose.tools import set_trace

class Derivator(object):

    def __init__(self, scalar_arr, len1, len2):
        self.arr = scalar_arr
        self.len1 = float(len1)
        self.len2 = float(len2)

        self.deriv1 = vcalc.cderivative(self.arr, 'X_DIR')
        self.deriv2 = vcalc.cderivative(self.arr, 'Y_DIR')

        self.perp_deriv1 =  self.deriv2
        self.perp_deriv2 = -self.deriv1

        deriv12 = vcalc.cderivative(self.deriv1, 'Y_DIR')
        deriv11 = vcalc.cderivative(self.arr, 'X_DIR', order=2)
        deriv22 = vcalc.cderivative(self.arr, 'Y_DIR', order=2)

        self.deriv12_interp = Interp2DPeriodic(len1, len2, deriv12)
        self.deriv11_interp = Interp2DPeriodic(len1, len2, deriv11)
        self.deriv22_interp = Interp2DPeriodic(len1, len2, deriv22)

class NullCell(object):

    def __init__(self, deriv, loc, size):
        self.deriv = deriv
        if size not in (1,2):
            raise ValueError('argument size must be in (1,2)')
        self.loc = loc
        self.size = size
        self.bounds = self.deriv.perp_deriv1.shape

    def find_null(self):
        nx,ny = self.bounds
        x,y = self.loc
        pt1 = x, y
        pt2 = x, (y+self.size) % ny
        pt3 = (x+self.size) % nx, (y+self.size) % ny
        pt4 = (x+self.size) % nx, y

        zf = []
        for pair in ((pt1, pt2), (pt2, pt3), (pt4, pt3), (pt1, pt4)):
            l0, l1 = pair
            x0, x1 = self.deriv.perp_deriv1[l0], self.deriv.perp_deriv1[l1]
            y0, y1 = self.deriv.perp_deriv2[l0], self.deriv.perp_deriv2[l1]
            zero_x = zero_loc(x0, x1)
            zero_y = zero_loc(y0, y1)
            if 0.0 <= zero_x <= 1.0:
                zero_l = zero_x
            elif 0.0 <= zero_y <= 1.0:
                zero_l = zero_y
            else:
                raise RuntimeError("unable to find zero crossing")
            zf.append(zero_l)

        #   Y   0 ---------------> 1
        #  X
        #  0    1 --------A------- 2  
        #  |    |         |        |
        #  |    |         |        |
        #  |    D - - - - + - - - -B
        #  |    |         |        |
        #  v    |         |        |
        #  1    4 --------C------- 3

        # A = (0.0, zf[0])
        # B = (zf[1], 1.0)
        # C = (1.0, zf[2])
        # D = (zf[3], 0.0)

        # We solve the linear system

        # V_1 = A - alpha (A - C)
        # V_2 = D + beta  (B - D)

        # for alpha and beta, such that

        # V_1(alpha) == V_2(beta).

        # Compute the x and y differences in opposing points.
        ACx = -1.0
        BDx = zf[1] - zf[3]
        ACy = zf[0] - zf[2]
        BDy = 1.0
        ADx = -zf[3]
        ADy =  zf[0]

        # Solve M (alpha beta).T = (ADx ADy).T
        # where M == [[ACx BDx] [ACy BDy]]

        # Hardcode M_inverse
        det = ACx * BDy - ACy * BDx
        M_inv = np.array([[BDy, -BDx], [-ACy, ACx]]) / det

        alpha, _ = np.dot(M_inv, [ADx, ADy])

        # Substitute the solved alpha into the eqn for V_1.
        # Yields the intersection point (x_soln, y_soln)
        x_soln = 0.0 - alpha * ACx
        y_soln = zf[0] - alpha * ACy

        # put the solution point in global coordinates.
        return (x + self.size * x_soln, y + self.size * y_soln)

    def is_saddle(self):
        x0, y0 = self.find_null()
        evals, evecs = eigsystem(self.deriv, x0, y0)
        return null_is_saddle(evals)

def zero_loc(y0, y1):
    return -y0 / (y1 - y0)

def flip_lr(p1, p2, p3, p4):
    return p1 * p2 <= 0.0e0 and p3 * p4 <= 0.0e0

def flip_tb(p1, p2, p3, p4):
    return p1 * p3 <= 0.0e0 and p2 * p4 <= 0.0e0

def is_null(x1, x2, x3, x4, y1, y2, y3, y4):
    return (flip_lr(x1, x2, x3, x4) and flip_tb(y1, y2, y3, y4) or
        flip_lr(y1, y2, y3, y4) and flip_tb(x1, x2, x3, x4))

def find_and_cull_cells(deriv):
    null_step2 = find_null_cells(deriv, step=2)
    null_step1 = find_null_cells(deriv, step=1)
    all_nulls = null_step2 + null_step1
    return remove_overlaps(all_nulls)

def find_null_cells(deriv, step=1):
    null_cells = []
    nx, ny = deriv.perp_deriv1.shape
    assert (nx, ny) == deriv.perp_deriv2.shape
    for i in range(nx):
        for j in range(ny):

            # # load up the x's and y's
            # #
            # #  y ---->
            # #  1 --- 2  x
            # #  |     |  |
            # #  |     |  |
            # #  3 --- 4  V

            istep = (i+step) % nx
            jstep = (j+step) % ny

            xb1 = deriv.perp_deriv1[i, j]
            xb2 = deriv.perp_deriv1[i, jstep]
            xb3 = deriv.perp_deriv1[istep, j]
            xb4 = deriv.perp_deriv1[istep, jstep]

            yb1 = deriv.perp_deriv2[i, j]
            yb2 = deriv.perp_deriv2[i, jstep]
            yb3 = deriv.perp_deriv2[istep, j]
            yb4 = deriv.perp_deriv2[istep, jstep]

            if is_null(xb1, xb2, xb3, xb4, yb1, yb2, yb3, yb4):
                null_cells.append(NullCell(deriv, (i,j), step))

    return null_cells

def remove_overlaps(cells):
    nx, ny = cells[0].bounds
    loc2cells = defaultdict(list)
    for cell in cells:
        loc2cells[cell.loc].append(cell)

    cell2cluster = defaultdict(set)
    for cell in cells:
        cell2cluster[cell].add(cell)

    offsets= ((-1,-1), (-1,0), (-1,1),
              ( 0,-1), (0, 0), ( 0,1),
              ( 1,-1), (1, 0), ( 1,1))

    for cell in cells:
        # get the cell's neighbors
        neighbors = []
        for offset in offsets:
            oi, oj = offset
            i, j = cell.loc
            nbr_loc = (oi+i) % nx, (oj+j) % ny
            if nbr_loc in loc2cells:
                neighbors.extend(loc2cells[nbr_loc])
        neighbors.remove(cell)
        # update the cell2cluster mapping
        cluster = cell2cluster[cell]
        for nbr in neighbors:
            cluster.update(cell2cluster[nbr])
            cell2cluster[nbr] = cluster

    clusters = dict((id(v), v) for v in cell2cluster.values()).values()

    reduced = []

    for idx, cluster in enumerate(clusters[:]):
        onesize = [cell for cell in cluster if cell.size == 1]
        assert len(onesize) in (0,1)
        if onesize:
            reduced.extend(onesize)
        else:
            sorted_cells = sorted(cluster, key=lambda x: x.loc)
            reduced.append(sorted_cells[0])

    return reduced

def eigsystem(deriv, x1, x2):
    'The psi_NN are interpolator objects to be evaluated at position (x1, x2).'

    from numpy.linalg import eig

    psi_12_val = deriv.deriv12_interp.eval(x1, x2)
    lin_matrix = np.array(
            [[psi_12_val,            deriv.deriv22_interp.eval(x1, x2)],
             [-deriv.deriv11_interp.eval(x1, x2), -psi_12_val         ]])
    evals, evecs = eig(lin_matrix)
    return evals, evecs

def null_is_saddle(evals):
    '''
    Inspects the given 2 eigenvalues, returns True iff evals are real and add
    to zero <==> null point is a saddle point.  Returns False otherwise.

    'Circulation' nulls have complex eigenvalues. Eigenvalues are purely
    imaginary when the circulation null is neither source or sink.

    '''
    if np.any(np.iscomplex(evals)):
        assert np.allclose(evals[0], evals[1].conjugate())
        return False
    if np.all(np.isreal(evals)) and np.allclose(evals[0], -evals[1]):
        return True
    raise ValueError("an interesting exception, evals==%s", list(evals))
