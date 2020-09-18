import numpy as np

def analytic_solution(cs):
    Ez = np.zeros((cs.nt, cs.nx))
    By = np.zeros((cs.nt, cs.nx))

    for n in range(0, cs.nt):
        t = n * cs.dt

        # k x E = w B
        Ez[n, 0] = np.sin(cs.w * t)
        By[n, 0] = (-1 / cs.imp0) * np.sin(cs.w * t)

        # Move wave forward one step
        Ez[n, 1:] = Ez[n - 1, :-1]
        By[n, 1:] = By[n - 1, :-1]

    return Ez, By
