from plotting import *
from analytic import *
from explicit import *
from implicit import *


class Constants:
    def __init__(self, x_res, num_l, cfl):
        self.c = 299792458   # m/s
        self.mu = 1.25663706e-6
        self.eps = 8.85418782e-12

        self.f = 14.0e6          # frequency of 14 MHz
        self.l = self.c / self.f           # wavelength
        self.w = 2 * np.pi * self.f   # omega
        self.imp0 = 377.0        # impedance of free space

        self.x_res = x_res
        self.num_l = num_l
        self.cfl = cfl

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


def run_sim(cs, file1, file2, plot=False, save=False):
    Ez_ass, By_ass = analytic_solution(cs)
    Ez_yee, By_yee = explicit_yee(cs)
    Ez_imp, By_imp = implicit_solution(cs)

    # ez_exp_err, ez_imp_err, ez_diff_err = error(Ez_ass[-1, :], Ez_yee[-1, :], Ez_imp[-1, :])
    # by_exp_err, by_imp_err, by_diff_err = error(By_ass[-1, :], By_yee[-1, :], By_imp[-1, :])
    #
    # print(f'\t\t\t\t   Ez\t\t   By')
    # print('-------------------------------------')
    # print(f'Explicit Error: {ez_exp_err:.3E}\t{by_exp_err:.3E}')
    # print(f'Implicit Error: {ez_imp_err:.3E}\t{by_imp_err:.3E}')
    # print(f'Imp-Exp  Error: {ez_diff_err:.3E}\t{by_diff_err:.3E}\n')

    if plot:
        multi_animated_lines([Ez_ass, Ez_yee, Ez_imp], cs.nt, file1, save=save)
        multi_animated_lines([By_ass, By_yee, By_imp], cs.nt, file2, save=save)


def main():
    exp_consts = Constants(16, 2, 1.0)
    imp_consts = Constants(32, 2, 2.0)
    save = False
    # Filenames to save animated plots
    e_filename = 'Ez_overlayed.mp4'
    b_filename = 'By_overlayed.mp4'

    Ez_ass, By_ass = analytic_solution(exp_consts)
    Ez_exp, By_exp = explicit_yee(exp_consts)
    Ez_imp, By_imp = implicit_solution(imp_consts)

    low_nt = min(exp_consts.nt, imp_consts.nt)

    multi_animated_lines([Ez_ass, Ez_exp, Ez_imp], low_nt, e_filename, save=save)
    multi_animated_lines([By_ass, By_exp, By_imp], low_nt, b_filename, save=save)


if __name__ == '__main__':
    main()
