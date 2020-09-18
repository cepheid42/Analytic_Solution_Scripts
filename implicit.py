import numpy as np

# Tridagonal Solver / Thomas Algorithm
def tridiagonal(d, cs):
    coeff = (1.0 / (8 * cs.eps * cs.mu)) * (cs.dt / cs.dx) ** 2

    nn = len(d)

    aa = np.full(nn - 1, -coeff)
    bb = np.full(nn, 0.5 + 2 * coeff)

    A = np.diag(bb) + np.diag(aa, 1) + np.diag(aa, -1)

    return np.linalg.solve(A, d)

    # for j in range(1, nn):
    #     mc = aa[j - 1] / bb[j - 1]
    #     bb[j] = bb[j] - mc * cc[j - 1]
    #     dd[j] = dd[j] - mc * dd[j - 1]
    #
    # x = bb
    # x[-1] = dd[-1] / bb[-1]
    # for j in range(nn - 2, -1, -1):
    #     x[j] = (dd[j] - cc[j] * x[j + 1]) / bb[j]
    #
    # return np.asarray(x)

def bp():
    pass


def implicit_solution(cs):
    Ez = np.zeros(cs.nx)
    By = np.zeros(cs.nx)

    Ez_out = [np.copy(Ez)]
    By_out = [np.copy(By)]

    c1 = (cs.dt / (2 * cs.eps * cs.dx))      # e-field coefficient
    c2 = (cs.dt / (2 * cs.mu * cs.dx))       # b-field coefficient

    for n in range(cs.nt - 1):
        t = n * cs.dt

        # Source
        Ez[0] = np.sin(cs.w * t)

        for x in range(cs.nx):
            By[x] = By[x] + c2 * (Ez[x + 1] - Ez[x - 1])

        bp()

        for x in range(1, cs.nx - 1):
            Ez[x] = Ez[x] + c1 * (By[x + 1] - By[x - 1])

        Ez_out.append(np.copy(Ez))
        By_out.append(np.copy(By))

    return np.asarray(Ez_out), np.asarray(By_out)
