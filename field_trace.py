import numpy as np
from nose.tools import set_trace

from collections import defaultdict

class NullCell(object):

    def __init__(self, ax, ay, loc, size):
        if size not in (1,2):
            raise ValueError('argument size must be in (1,2)')
        assert ax.shape == ay.shape
        self.ax, self.ay = ax, ay
        self.loc = loc
        self.size = size
        self.bounds = ax.shape
        # self.null_loc = self.find_null()

    def covered(self):
        if self.size == 1:
            return set((self.loc,))
        if self.size == 2:
            offsets = ((0,0), (0,1), (1,0), (1,1))
            covered = set()
            i,j = self.loc
            nx, ny = self.bounds
            for offset in offsets:
                oi, oj = offset
                new_loc = ((oi+i) % nx, (oj+j) % ny)
                covered.add(new_loc)
            return covered

    def find_null(self):
        nx,ny = self.bounds
        x,y = self.loc
        pt1 = self.loc
        pt2 = x, (y+self.size) % ny
        pt3 = (x+self.size) % nx, (y+self.size) % ny
        pt4 = (x+self.size) % nx, y

        zf = []
        for pair in ((pt1, pt2), (pt2, pt3), (pt4, pt3), (pt1, pt4)):
            l0, l1 = pair
            x0, x1 = self.ax[l0], self.ax[l1]
            y0, y1 = self.ay[l0], self.ay[l1]
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

def slope_intercept(x1, y1, x2, y2):
    slope = (y2 - y1) / (x2 - x1)
    intcpt = y1 - slope * x1
    return (slope, intcpt)

def zero_loc(y0, y1):
    return -y0 / (y1 - y0)

def flip_lr(p1, p2, p3, p4):
    return p1 * p2 <= 0.0e0 and p3 * p4 <= 0.0e0

def flip_tb(p1, p2, p3, p4):
    return p1 * p3 <= 0.0e0 and p2 * p4 <= 0.0e0

def is_null(x1, x2, x3, x4, y1, y2, y3, y4):
    return (flip_lr(x1, x2, x3, x4) and flip_tb(y1, y2, y3, y4) or
        flip_lr(y1, y2, y3, y4) and flip_tb(x1, x2, x3, x4))

def find_and_cull_cells(ax, ay):
    null_step2 = find_null_cells(ax, ay, step=2)
    null_step1 = find_null_cells(ax, ay, step=1)
    all_nulls = null_step2 + null_step1
    return remove_overlaps(all_nulls)

def find_null_cells(ax, ay, step=1):
    null_cells = []
    nx, ny = ax.shape
    assert (nx, ny) == ay.shape
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

            xb1 = ax[i, j]
            xb2 = ax[i, jstep]
            xb3 = ax[istep, j]
            xb4 = ax[istep, jstep]

            yb1 = ay[i, j]
            yb2 = ay[i, jstep]
            yb3 = ay[istep, j]
            yb4 = ay[istep, jstep]

            if is_null(xb1, xb2, xb3, xb4, yb1, yb2, yb3, yb4):
                null_cells.append(NullCell(ax, ay, (i,j), step))

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

def _remove_overlaps(cells):
    '''
    given a list of cells, removes cells from the list such that no 2 cells in
    the returned list overlap.

    2 cells a & b overlap iff len(a.covered().intersection(b.covered())) > 0
    '''

    covered2cells = defaultdict(list)
    for cell in cells:
        for cov in cell.covered:
            covered2cells[cov].append(cell)

    for cov in covered2cells:
        cells = covered2cells[cov]
        # Preference for cells with size == 1, since these are more specific.
        onesize = [cell for cell in cells if cell.size == 1]
        assert len(onesize) in (0,1)
        if onesize:
            # remove all other cells, keeping just this cell with size == 1.
            covered2cells[cov] = onesize
        else:
            # sort the size-2 cells based on cell.loc and pick the first in the list.
            sorted_cells = sorted(cells, key=lambda x: x.loc)
            covered2cells[cov] = sorted_cells[0]

    uniq_cells = set(covered2cells.values())
    return list(uniq_cells)

def eigsystem(psi_12, psi_22, psi_11, x1, x2):
    'The psi_NN are interpolator objects to be evaluated at position (x1, x2).'

    from numpy.linalg import eig

    psi_12_val = psi_12.eval(x1, x2)
    lin_matrix = np.array(
            [[psi_12_val,            psi_22.eval(x1, x2)],
             [-psi_11.eval(x1, x2), -psi_12_val         ]])
    evals, evecs = eig(lin_matrix)
    return evals, evecs

def null_is_saddle(evals):
    if np.any(np.iscomplex(evals)):
        assert np.allclose(evals[0], evals[1].conjugate())
        return False
    if np.all(np.isreal(evals)) and np.allclose(evals[0], -evals[1]):
        return True
    raise ValueError("an interesting exception, evals==%s", list(evals))
