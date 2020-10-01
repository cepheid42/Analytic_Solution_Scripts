import numpy as np

# Tridagonal Solver / Thomas Algorithm
def tridiagonal(d, cs):
    coeff = (1.0 / (8 * cs.eps * cs.mu)) * (cs.dt / cs.dx) ** 2

    aa = np.full(cs.nx - 1, -coeff)
    bb = np.full(cs.nx, 0.5 + 2 * coeff)

    A = np.diag(bb) + np.diag(aa, 1) + np.diag(aa, -1)

    return np.linalg.solve(A, d)

def bp():
    pass


def implicit_solution(cs):
    Ez = np.zeros((cs.nt, cs.nx))
    By = np.zeros((cs.nt, cs.nx))

    # c1 = (cs.dt / (4 * cs.eps * cs.dx))      # e-field coefficient
    # c2 = (cs.dt / (4 * cs.mu * cs.dx))       # b-field coefficient
    #
    # for n in range(1, cs.nt):
    #     t = n * cs.dt
    #     ez_rhs = np.zeros(cs.nx)
    #     ez_rhs[0] = Ez[n - 1, 0]
    #     ez_rhs[-1] = Ez[n - 1, -1]
    #     # Implicit e-field step
    #     for x in range(1, cs.nx - 1):
    #         ez_rhs[x] = Ez[n - 1, x] - c1 * (By[n - 1, x + 1] - By[n - 1, x - 1])
    #         if ez_rhs[x] > 1.0:
    #             bp()
    #
    #     # Tridiagonal solution
    #     ez = tridiagonal(ez_rhs, cs)
    #
    #     # Explicit Electric field step
    #     for x in range(cs.nx):
    #         Ez[n, x] = ez[x] - Ez[n - 1, x]
    #         if Ez[n, x] > 1.0:
    #             bp()
    #
    #     # Explicit Magnetic field step
    #     for x in range(1, cs.nx - 1):
    #         By[n, x] = By[n - 1, x] + c2 * (ez[x + 1] - ez[x - 1])
    #
    #         if By[n, x] > 1.0:
    #             bp()
    #
    #     # Source
    #     Ez[n, 0] = np.sin(cs.w * t)

    return Ez, By
