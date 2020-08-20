from plotting import *

c = 299792458   # m/s
mu = 1.25663706e-6
eps = 8.85418782e-12

f = 14.0e6      # 14 MHz
l = c / f
w = 2 * np.pi * f

x_res = 128
num_wavelengths = 2

nx = x_res * num_wavelengths + 1
dx = (num_wavelengths * l) / (nx - 1)

cfl = 1.0

dt = cfl * dx / c
nt = int(num_wavelengths * l / (c * dt)) + 1


def analytic_solution():
    waves = np.zeros((nt, nx))

    for n in range(1, nt):
        t = n * dt
        # Source
        waves[n, 0] = np.sin(w * t)
        # Move wave forward one step
        waves[n, 1:] = waves[n - 1, :-1]

    return waves


def explicit_yee():
    Ez = np.zeros((nt, nx))
    By = np.zeros((nt, nx))

    c1 = dt / (eps * dx)
    c2 = dt / (mu * dx)

    for n in range(1, nt):
        t = n * dt
        # Magnetic Update
        By[n, :-1] = By[n - 1, :-1] + c2 * (Ez[n - 1, 1:] - Ez[n - 1, :-1])
        # Electric Update
        Ez[n, 1:] = Ez[n - 1, 1:] + c1 * (By[n, 1:] - By[n, :-1])

        # Source
        Ez[n, 0] = np.sin(w * t)
    return Ez, By


def implicit_solution():
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

    Ez = np.zeros((nt, nx))
    By = np.zeros((nt, nx))

    c1 = (dt / dx) * (1 / (4 * eps))
    c2 = (dt / dx) * (1 / (4 * mu))

    for n in range(1, nt):
        t = n * dt

        # implicit Ez_half
        ez_rhs = Ez[n - 1, 1:-1] + c1 * (By[n - 1, 2:] - By[n - 1, :-2])
        ez_rhs = np.insert(ez_rhs, 0, Ez[n - 1, 0])
        ez_rhs = np.insert(ez_rhs, -1, Ez[n - 1, -1])
        ez = tridiagonal(ez_rhs)

        # Explicit Ez
        Ez[n, :] = ez - Ez[n - 1, :]

        # Explicit By
        By[n, 0] = By[n - 1, 0]
        By[n, -1] = By[n - 1, -1]
        By[n, 1:-1] = By[n - 1, 1:-1] + c2 * (ez[2:] - ez[:-2])

        # Source
        Ez[n, 0] = np.sin(w * t)

    return Ez, By


if __name__ == '__main__':
    analytic = analytic_solution()
    # animated_line(analytic, nt, "Analytic Solution")
    # waterfall(analytic, nx, nt)

    Ez_yee, By_yee = explicit_yee()
    # animated_line(Ez_yee, nt, "Explicit Yee Solution")
    # animated_line(By_yee, nt, "Explicit Yee Solution")

    Ez_imp, By_imp = implicit_solution()
    # animated_line(Ez_exp, nt, "Implicit Ez")
    # animated_line(By_exp, nt, "Implicit By")

    file1 = f'Ez_multi_{x_res}_{num_wavelengths}.mp4'
    # file2 = f'By_multi_{x_res}_{num_wavelengths}.mp4'
    multi_animated_lines([analytic, Ez_yee, Ez_imp], nt, file1) #, save=True)
    # multi_animated_lines([analytic, By_yee, By_imp], nt, file2, save=True)
