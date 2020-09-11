from plotting import *
from analytic import *
from explicit import *
from implicit import *

class Constants:
    def __init__(self):
        self.c = 299792458   # m/s
        self.mu = 1.25663706e-6
        self.eps = 8.85418782e-12

        self.f = 14.0e6          # frequency of 14 MHz
        self.l = self.c / self.f           # wavelength
        self.w = 2 * np.pi * self.f   # omega
        self.imp0 = 377.0        # impedance of free space

        self.x_res = 64
        self.num_l = 2
        self.cfl = 1.0

        self.nx = self.x_res * self.num_l + 1
        self.dx = (self.num_l * self.l) / (self.nx - 1)
        self.dt = self.cfl * self.dx / self.c
        self.nt = int(self.num_l * self.l / (self.c * self.dt)) + 1

def error(analytic, explicit, implicit):
    # Mean Squared Error
    exp_err = np.mean((explicit - analytic)**2)
    imp_err = np.mean((implicit - analytic)**2)
    diff_err = np.mean((implicit - explicit)**2)

    return exp_err, imp_err, diff_err

def run_sim(nt, dt, nx, dx, file1, file2, plot=False, save=False):
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

def main():
    consts = Constants()

    # Filenames to save animated plots
    e_filename = 'Ez_overlayed.mp4'
    b_filename = 'By_overlayed.mp4'

    # Ey_ass, Ez_ass, By_ass, Bz_ass = analytic_solution(consts)
    # animated_line(Ez_ass, By_ass, consts, 'z', "Ez-By Analytic")
    # animated_line(Ey_ass, Bz_ass, consts, 'y', "Ey-Bz Analytic")

    Ey_exp, Ez_exp, By_exp, Bz_exp = explicit_yee(consts)
    animated_line(Ez_exp, By_exp, consts, 'z', "Ez-By Explicit")
    animated_line(Ey_exp, By_exp, consts, 'y', "Ey-Bz Explicit")


if __name__ == '__main__':
    main()
