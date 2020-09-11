import numpy as np

def analytic_solution(cs):
    Ey = np.zeros((cs.nt, cs.nx), dtype=np.float64)
    Ez = np.zeros((cs.nt, cs.nx), dtype=np.float64)

    By = np.zeros((cs.nt, cs.nx), dtype=np.float64)
    Bz = np.zeros((cs.nt, cs.nx), dtype=np.float64)

    for n in range(0, cs.nt):
        t = n * cs.dt
        # k x E = w B
        # Ex = Bx = 0
        # kx * Ez = -w By
        # kx * Ey = w Bz
        Ey[n, 0] = np.sin(cs.w * t)
        Ez[n, 0] = np.sin(cs.w * t)

        By[n, 0] = (-1 / cs.imp0) * np.sin(cs.w * t)
        Bz[n, 0] = (1 / cs.imp0) * np.sin(cs.w * t)

        # Move wave forward one step
        Ey[n, 1:] = Ey[n - 1, :-1]
        Ez[n, 1:] = Ez[n - 1, :-1]

        By[n, 1:] = By[n - 1, :-1]
        Bz[n, 1:] = Bz[n - 1, :-1]

    return Ey, Ez, By, Bz
