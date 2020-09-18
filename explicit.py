import numpy as np


def explicit_yee(cs):
    Ez = np.zeros((cs.nt, cs.nx))
    By = np.zeros((cs.nt, cs.nx))

    # c1 = cs.dt / (eps * dx)
    # c2 = cs.dt / (mu * dx)

    for n in range(1, cs.nt):
        t = n * cs.dt

        # Magnetic Update
        By[n, :-1] = By[n - 1, :-1] + (Ez[n - 1, 1:] - Ez[n - 1, :-1]) / cs.imp0

        # Electric Update
        Ez[n, 1:] = Ez[n - 1, 1:] + (By[n, 1:] - By[n, :-1]) * cs.imp0

        # Source
        Ez[n, 0] = np.sin(cs.w * t)

    return Ez, By
