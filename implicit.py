import numpy as np


# Tridagonal Solver / Thomas Algorithm
def tridiagonal(d, cs):
    coeff = 1.0 / (8 * cs.eps * cs.mu) * (cs.dt / cs.dx) ** 2
    dd = np.copy(d)
    nn = len(d)

    aa = np.full(nn - 1, -coeff)
    cc = np.full(nn - 1, -coeff)
    bb = np.full(nn, 0.5 + 2 * coeff)

    for j in range(1, nn):
        mc = aa[j - 1] / bb[j - 1]
        bb[j] = bb[j] - mc * cc[j - 1]
        dd[j] = dd[j] - mc * dd[j - 1]

    x = bb
    x[-1] = dd[-1] / bb[-1]
    for j in range(nn - 2, -1, -1):
        x[j] = (dd[j] - cc[j] * x[j + 1]) / bb[j]

    return np.asarray(x)


def implicit_solution(cs):
    Ey = np.zeros((cs.nt, cs.nx))
    Ez = np.zeros((cs.nt, cs.nx))

    By = np.zeros((cs.nt, cs.nx))
    Bz = np.zeros((cs.nt, cs.nx))

    c1 = (cs.dt / (4 * cs.eps * cs.dx))      # e-field coefficient
    c2 = (cs.dt / (4 * cs.mu * cs.dx))       # b-field coefficient

    for n in range(0, cs.nt):
        t = n * cs.dt

        # implicit step
        ez_rhs = Ez[n - 1, 1:-1] + c1 * (By[n - 1, 2:] - By[n - 1, :-2])

        # These are needed to fix the boundaries, just inserts the left and right edges
        # from the last time step into the array
        ez_rhs = np.insert(ez_rhs, 0, Ez[n - 1, 0])
        ez_rhs = np.insert(ez_rhs, -1, Ez[n - 1, -1])

        ez = tridiagonal(ez_rhs)

        # Explicit Ez
        Ez[n, :] = ez - Ez[n - 1, :]

        # Explicit By
        By[n, 1:-1] = By[n - 1, 1:-1] + c2 * (ez[2:] - ez[:-2])
        # carry left and right edges forward one time step
        By[n, 0] = By[n - 1, 0]
        By[n, -1] = By[n - 1, -1]

        # Source
        Ez[n, 0] = np.sin(cs.w * t)

    return Ez, By