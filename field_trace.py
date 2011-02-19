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

def find_null(deriv, i, j):

    a1, b1, c1, d1 = _get_abcd(deriv.perp_deriv1, i, j)
    a2, b2, c2, d2 = _get_abcd(deriv.perp_deriv2, i, j)

    if 0:
        # First version:

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

    if 1:
        # Second version:
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


class NullCell(object):

    def __init__(self, deriv, loc, size):
        self.deriv = deriv
        if size not in (1,2):
            raise ValueError('argument size must be in (1,2)')
        self.loc = loc
        self.size = size
        self.bounds = self.deriv.perp_deriv1.shape
        self.x0, self.y0 = find_null(self.deriv, *self.loc)
        self.evals, self.evecs = eigsystem(self.deriv, self.x0, self.y0)

    def _find_null(self):
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
                import pdb; pdb.set_trace()
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
        return null_is_saddle(self.evals)

    def incoming_disp(self):
        if self.evals[0] < 0.0:
            return self.evecs[:,0]
        elif self.evals[1] < 0.0:
            return self.evecs[:,1]

    def outgoing_disp(self):
        if self.evals[0] > 0.0:
            return self.evecs[:,0]
        elif self.evals[1] > 0.0:
            return self.evecs[:,1]

    def outgoing_start(self, scale):
        if self.evals[0] > 0.0:
            evec = self.evecs[:,0]
        elif self.evals[1] > 0.0:
            evec = self.evecs[:,1]
        return ((self.x0 + scale * evec[0],
                 self.y0 + scale * evec[1]),
                (self.x0 - scale * evec[0],
                 self.y0 - scale * evec[1]))

def zero_loc(y0, y1):
    return -y0 / (y1 - y0)

def flip_lr(p1, p2, p3, p4):
    return p1 * p2 <= 0.0e0 and p3 * p4 <= 0.0e0

def flip_tb(p1, p2, p3, p4):
    return p1 * p3 <= 0.0e0 and p2 * p4 <= 0.0e0

def is_null_flip_test(x1, x2, x3, x4, y1, y2, y3, y4):
    return (flip_lr(x1, x2, x3, x4) and flip_tb(y1, y2, y3, y4) or
        flip_lr(y1, y2, y3, y4) and flip_tb(x1, x2, x3, x4))

def same_sign_or_zero(*args):
    for arg in args:
        if arg:
            test_sgn = -1 if arg < 0 else 1
    sgns = np.sign(args)
    return np.all((sgns == test_sgn) | (sgns == 0))

def same_sign(*args):
    sgns = np.sign(args)
    return np.all(sgns == sgns[0]) and sgns[0] in (1, -1)

def is_null(xs, ys):
    if same_sign(*xs) or same_sign(*ys):
        return False
    return True

def find_and_cull_cells(deriv):
    # null_step2 = find_null_cells(deriv, step=2)
    null_step1 = find_null_cells(deriv, step=1)
    return null_step1
    # all_nulls = null_step2 + null_step1
    # return remove_overlaps(all_nulls)

def _get_corner_vals(f, i, j):

    # #  y --------->
    # #  
    # #  f00 ---- f01  x
    # #  |         |   |
    # #  |         |   |
    # #  f10 ---- f11  V

    nx, ny = f.shape
    step = 1
    istep = (i+step) % nx
    jstep = (j+step) % ny

    x00 = f[i, j]
    x01 = f[i, jstep]
    x10 = f[istep, j]
    x11 = f[istep, jstep]

    return (x00, x01, x10, x11)

def _get_abcd(f, i, j):
    f00, f01, f10, f11 = _get_corner_vals(f, i, j)

    a = f00
    b = f10 - f00
    c = f01 - f00
    d = f11 - f10 - f01 + f00

    return (a, b, c, d)

def find_null_cells(deriv, step=1):
    step = 1
    null_cells = []
    nx, ny = deriv.perp_deriv1.shape
    assert (nx, ny) == deriv.perp_deriv2.shape
    for i in range(nx):
        for j in range(ny):
            if is_null(_get_corner_vals(deriv.perp_deriv1, i, j),
                       _get_corner_vals(deriv.perp_deriv2, i, j)):
                try:
                    find_null(deriv, i, j)
                except ValueError:
                    pass
                else:
                    null_cells.append(NullCell(deriv, (i,j), step))

    return null_cells

'''
def _find_null_cells(deriv, step=1):
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
'''

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
        if len(onesize) not in (0,1):
            # FIXME: this should be an error... why??
            import pdb; pdb.set_trace()
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

def make_skeleton(arr, deriv, saddles, ncomb, start_offset, hit_radius, timeout=500):
    from wrap_trace_integrator import TraceIntegrator

    scale = 1. / max(deriv.deriv1.max(), deriv.deriv2.max())

    ti = TraceIntegrator(ntracers=ncomb, scale=scale,
            x0max=deriv.len1, x1max=deriv.len2,
            v0=deriv.perp_deriv1_interp, v1=deriv.perp_deriv2_interp,
            eps_abs=1.0e-6, eps_rel=0.0)

    def within_hit_rad(comb_positions, saddle_positions):
        comb_center = comb_positions.sum(axis=0) / comb_positions.shape[0]
        comb_radius2 = np.max(np.sum((comb_positions - comb_center)**2, axis=1))
        distances2 = np.sum((saddle_positions - comb_center)**2, axis=1)
        if np.any(distances2 <= comb_radius2):
            return saddle_positions[distances2 <= comb_radius2]
        else:
            return []
        
    comb_positions = get_comb(saddles[0], ncomb, start_offset, hit_radius)

    saddle_positions = np.array(
                            [[s.x0, s.y0] for s in saddles if s is not saddles[0]])

    dt = 1. / 2**3

    trace = []

    t_init = 0.0
    ctr = 0
    while True:
        ctr += 1
        trace.append(list(comb_positions.flatten()))
        ti.evolve(t=t_init, t1=t_init+dt, h=1.0e-6, y=comb_positions.ravel())
        t_init += dt
        clip_positions(comb_positions, deriv.len1, deriv.len2)
        if t_init >= timeout:
            return (False, trace)
        hits = within_hit_rad(comb_positions, saddle_positions)
        if np.any(hits):
            return (True, trace)

def trace_from_saddle(arr, deriv, start_saddle, start_offset, outgoing=True, timeout=500):
    from wrap_trace_integrator import TraceIntegrator

    scale = 10. / max(deriv.deriv1.max(), deriv.deriv2.max())
    if not outgoing:
        scale *= -1

    ti = TraceIntegrator(ntracers=1, scale=scale,
            x0max=deriv.len1, x1max=deriv.len2,
            v0=deriv.perp_deriv1_interp, v1=deriv.perp_deriv2_interp,
            eps_abs=1.0e-6, eps_rel=0.0)

    offset = np.array((start_saddle.x0, start_saddle.y0))
    if outgoing:
        evec = start_saddle.outgoing_disp()
    else:
        evec = start_saddle.incoming_disp()
    position = np.array([offset + start_offset * evec])

    dt = 1. / 2**1

    trace = []

    t_init = 0.0

    while True:
        trace.append(list(position.flatten()))
        ti.evolve(t=t_init, t1=t_init+dt, h=1.0e-6, y=position.ravel())
        t_init += dt
        clip_positions(position, deriv.len1, deriv.len2)
        if t_init >= timeout:
            return trace

def get_comb(nullcell, npts, start_offset, scale):
    offset = np.array((nullcell.x0, nullcell.y0))
    evec = nullcell.outgoing_disp()
    perp_evec = np.array([-evec[1], evec[0]])
    comb = np.array([(offset + start_offset * evec + ss * perp_evec) for ss in np.linspace(-scale, scale, npts, endpoint=True)])
    return comb

def clip_positions(posns, x0max, x1max):
    from wrap_trace_integrator import clip
    for pos in posns:
        pos[0] = clip(pos[0], x0max)
        pos[1] = clip(pos[1], x1max)

def neighbors(i, j, nx, ny):
    return (
            ((i-1)%nx, (j-1)%ny),
            ((i+1)%nx, (j-1)%ny),
            ((i+1)%nx, (j+1)%ny),
            ((i-1)%nx, (j+1)%ny),

            ((i-1)%nx, j),
            (i, (j-1)%ny),
            ((i+1)%nx, j),
            (i, (j+1)%ny)
            )

def level_sets(arr, interp, nulls):
    null2marked = {}
    for idx, null in enumerate(nulls):
        level_val = interp.eval(null.x0, null.y0)
        null2marked[null] = level_set(arr, level_val, (null.x0, null.y0))
    return null2marked

def marked_to_mask(shape, marked):
    mask = np.zeros(shape, dtype=np.bool_)
    for m in marked:
        for (i,j) in m:
            mask[i,j] = True
    return mask

def level_set(arr, level_val, position):
    from collections import deque
    nx, ny = arr.shape
    larr = arr - level_val # leveled array
    i0, j0 = int(position[0]) % nx, int(position[1]) % ny
    cvs = _get_corner_vals(larr, i0, j0)
    # FIXME: this should be the case, is it an interpolation issue?
    # if same_sign_or_zero(*cvs):
        # import pdb; pdb.set_trace()
    marked = set([(i0, j0)])
    front = deque([(i0, j0)])
    while front:
        i0, j0 = front.pop()
        for (i, j) in neighbors(i0, j0, nx, ny):
            if (i,j) in marked:
                continue
            cvs = _get_corner_vals(larr, i, j)
            if not same_sign_or_zero(*cvs):
                marked.add((i, j))
                front.appendleft((i, j))

    return list(marked)

def flood_fill(arr, i, j, border_color, fill_color):
    def neighbors(i, j, nx, ny):
        return (((i+1)%nx, j),
                (i, (j+1)%ny),
                ((i-1)%nx, j),
                (i, (j-1)%ny))

    nx, ny = arr.shape
    i %= nx; j %= ny
    if arr[i, j] == border_color:
        return
    edge = [(i,j)]
    arr[i,j] = fill_color
    while edge:
        newedge = []
        for (i, j) in edge:
            for (s, t) in neighbors(i, j, nx, ny):
                if arr[s,t] not in (border_color, fill_color):
                    arr[s,t] = fill_color
                    newedge.append((s,t))
        edge = newedge
