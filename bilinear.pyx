cpdef double bilinear(double f00,
                     double f10,
                     double f01,
                     double f11,
                     double x,
                     double y):
    return (f00 * (1. - x) * (1. - y) +
            f10 * x        * (1. - y) +
            f01 * (1. - x) * y        +
            f11 * x        * y)
