import numpy as np

def explicit_yee(cs):
    Ey = np.zeros((cs.nt, cs.nx))
    Ez = np.zeros((cs.nt, cs.nx))

    By = np.zeros((cs.nt, cs.nx))
    Bz = np.zeros((cs.nt, cs.nx))

    # c1 = cs.dt / (eps * dx)
    # c2 = cs.dt / (mu * dx)

    for n in range(1, cs.nt):
        t = n * cs.dt

        # Magnetic Update
        By[n, :-1] = By[n - 1, :-1] + (Ez[n - 1, 1:] - Ez[n - 1, :-1]) / cs.imp0
        Bz[n, :-1] = Bz[n - 1, :-1] - (Ey[n - 1, 1:] - Ey[n - 1, :-1]) / cs.imp0

        # Electric Update
        Ey[n, 1:] = Ey[n - 1, 1:] - (Bz[n, 1:] - Bz[n, :-1]) * cs.imp0
        Ez[n, 1:] = Ez[n - 1, 1:] + (By[n, 1:] - By[n, :-1]) * cs.imp0

        # Source
        Ez[n, 0] = np.sin(cs.w * t)
        Ey[n, 0] = np.sin(cs.w * t)

    return Ey, Ez, By, Bz
