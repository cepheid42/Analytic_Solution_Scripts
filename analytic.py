from plotting import *

c = 299792458   # m/s
mu = 1.25663706e-6
eps = 8.85418782e-12

f = 14.0e6          # frequency of 14 MHz
l = c / f           # wavelength
w = 2 * np.pi * f   # omega
imp0 = 377.0        # impedance of free space


def analytic_solution(nt, dt, nx):
    Ez = np.zeros((nt, nx), dtype=np.float64)
    By = np.zeros((nt, nx), dtype=np.float64)

    for n in range(0, nt):
        t = n * dt
        # Source is a transverse wave
        Ez[n, 0] = np.sin(w * t)
        By[n, 0] = (-1 / imp0) * np.sin(w * t)

        # Move wave forward one step
        Ez[n, 1:] = Ez[n - 1, :-1]
        By[n, 1:] = By[n - 1, :-1]

    return Ez, By


def explicit_yee(nt, dt, nx):
    Ez = np.zeros((nt, nx), dtype=np.float64)
    By = np.zeros((nt, nx), dtype=np.float64)

    # c1 = dt / (eps * dx)
    # c2 = dt / (mu * dx)

    for n in range(1, nt):
        t = n * dt
        # Magnetic Update
        By[n, :-1] = By[n - 1, :-1] + (Ez[n - 1, 1:] - Ez[n - 1, :-1]) / imp0
        # By[n, :-1] = By[n - 1, :-1] + c2 * (Ez[n - 1, 1:] - Ez[n - 1, :-1])

        # Electric Update
        Ez[n, 1:] = Ez[n - 1, 1:] + (By[n, 1:] - By[n, :-1]) * imp0
        # Ez[n, 1:] = Ez[n - 1, 1:] + c1 * (By[n, 1:] - By[n, :-1])

        # Source
        Ez[n, 0] = np.sin(w * t)

    return Ez, By


def implicit_solution(nt, dt, nx, dx):
    # Tridagonal Solver / Thomas Algorithm
    def tridiagonal(d):
        coeff = 1.0 / (8 * eps * mu) * (dt / dx) ** 2
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

    Ez = np.zeros((nt, nx), dtype=np.float64)
    By = np.zeros((nt, nx), dtype=np.float64)

    c1 = (dt / (4 * eps * dx))      # e-field coefficient
    c2 = (dt / (4 * mu * dx))       # b-field coefficient

    for n in range(0, nt):
        t = n * dt

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
        Ez[n, 0] = np.sin(w * t)

    return Ez, By


def error(analytic, explicit, implicit):
    """
    :param analytic: Array of analytic solution values
    :param explicit: Array of explicit solution values
    :param implicit: Array of implicit solution values
    :return: tuple of mean squared errors
    """
    # Mean Squared Error
    exp_err = np.mean((explicit - analytic)**2)
    imp_err = np.mean((implicit - analytic)**2)
    diff_err = np.mean((implicit - explicit)**2)

    return exp_err, imp_err, diff_err


def run_sim(nt, dt, nx, dx, file1, file2, plot=False, save=False):
    """
    :param file1: filename to save electric field
    :param file2: filename to save magnetic field
    :param plot: Plots results if True
    :param save: Save animated plots if True
    """

    Ez_ass, By_ass = analytic_solution(nt, dt, nx)
    Ez_yee, By_yee = explicit_yee(nt, dt, nx)
    Ez_imp, By_imp = implicit_solution(nt, dt, nx, dx)

    ez_exp_err, ez_imp_err, ez_diff_err = error(Ez_ass[-1, :], Ez_yee[-1, :], Ez_imp[-1, :])
    by_exp_err, by_imp_err, by_diff_err = error(By_ass[-1, :], By_yee[-1, :], By_imp[-1, :])

    print(f'\t\t\t\t   Ez\t\t   By')
    print('-------------------------------------')
    print(f'Explicit Error: {ez_exp_err:.3E}\t{by_exp_err:.3E}')
    print(f'Implicit Error: {ez_imp_err:.3E}\t{by_imp_err:.3E}')
    print(f'Imp-Exp  Error: {ez_diff_err:.3E}\t{by_diff_err:.3E}\n')

    if plot:
        multi_animated_lines([Ez_ass, Ez_yee, Ez_imp], nt, file1, save=save)
        multi_animated_lines([By_ass, By_yee, By_imp], nt, file2, save=save)


if __name__ == '__main__':
    x_res = 64              # number of points per wavelength
    n_wavelengths = 4
    cfl = 1.0               # Courant Number

    nx = x_res * n_wavelengths + 1
    dx = (n_wavelengths * l) / (nx - 1)

    dt = cfl * dx / c
    nt = int(n_wavelengths * l / (c * dt)) + 1

    # Filenames to save animated plots
    e_filename = 'Ez_overlayed.mp4'
    b_filename = 'By_overlayed.mp4'

    # Runs all three simulations
    run_sim(nt, dt, nx, dx, e_filename, b_filename, plot=True, save=False)
