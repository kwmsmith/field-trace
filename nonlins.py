from kaw_analysis.vcalc import gradient, laplacian

SIGMA = 8.0

class nonlins(object):

    def __init__(self, alpha, psi, phi, den):
        self.phi = phi
        self.psi = psi
        self.den = den
        self.alpha = alpha
        self.by, self.bx = gradient(psi)
        self.by *= -1
        self.vy, self.vx = gradient(phi)
        self.vy *= -1
        self.J = laplacian(psi)
        self.vor = laplacian(phi)

    def grad_par_b(self, arr):
        arr_x, arr_y = gradient(arr)
        return self.bx * arr_x + self.by * arr_y

    def grad_par_v(self, arr):
        arr_x, arr_y = gradient(arr)
        return self.vx * arr_x + self.vy *arr_y

    def psi_nonlin(self):
        # return self.alpha * self.grad_par_b(self.den) - self.grad_par_b(self.phi)
        return self.grad_par_b(self.alpha * self.den - self.phi)

    def vor_nonlin(self):
        return self.grad_par_v(self.vor) - self.grad_par_b(self.J)

    def den_nonlin(self):
        return self.grad_par_v(self.den) - self.grad_par_b(self.J)

    def grad_par_j(self):
        return self.grad_par_b(self.J)

    def grad_par_n(self):
        return self.alpha * self.grad_par_b(self.den)

