import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from matplotlib import animation


def waterfall(data, nx, nt):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    x_range = np.arange(nx - 1)

    count = 0
    for d in data:
        q = count * 2
        ax.plot(x_range, d[1:] + q, 'k', lw=1, zorder=nt - q)
        ax.fill_between(x_range, d[1:] + q, facecolor='w', lw=0, zorder=nt - q - 1)
        count += 1

    plt.show()


def animated_line(data, nt, title, save=False):
    fig = plt.figure()
    plts = []  # get ready to populate this list the Line artists to be plotted
    for i in range(nt):
        p, = plt.plot(data[i, :], 'k')  # this is how you'd plot a single line...
        plt.title(title)
        plts.append([p])  # ... but save the line artist for the animation
    ani = animation.ArtistAnimation(fig, plts, interval=50)  # run the animation
    if save:
        ani.save('wave.mp4')  # optionally save it to a file

    plt.show()


def tranverse_wave_plot(E, B, cs, e_dir, title, save=False):
    xs = np.linspace(0, cs.nx * cs.dx, E.shape[1])
    b_dir = 'By' if e_dir == 'z' else 'Bz'

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.view_init(elev=20., azim=-60)

    plts = []  # get ready to populate this list the Line artists to be plotted
    for i in range(cs.nt):
        if e_dir == 'z':
            ys1 = np.zeros(E.shape[1])
            zs1 = E[i, :]
            ys2 = B[i, :]
            zs2 = ys1
        else:
            ys1 = E[i, :]
            zs1 = np.zeros(E.shape[1])
            ys2 = zs1
            zs2 = B[i, :]

        pE, = ax.plot(xs=xs, ys=ys1, zs=zs1, zdir='z', color='r')  # this is how you'd plot a single line...
        pB, = ax.plot(xs=xs, ys=ys2, zs=zs2, zdir='z', color='b')
        plts.append([pE, pB])  # ... but save the line artist for the animation

    ani = animation.ArtistAnimation(fig, plts, interval=50, repeat=False)  # run the animation
    if save:
        ani.save(title + '.mp4')  # optionally save it to a file
    plt.legend(['E' + e_dir, b_dir])
    plt.show()


def multi_animated_lines(data, nt, title, save=False):
    fig = plt.figure()
    plts = []
    for i in range(nt):
        p1, = plt.plot(data[0][i, :], 'r')
        p2, = plt.plot(data[1][i, :], 'b--')
        p3, = plt.plot(data[2][i, :], 'g')
        plts.append([p1, p2, p3])  # ... but save the line artist for the animation
        # plts.append([p1, p2])
    ani = animation.ArtistAnimation(fig, plts, interval=50, repeat=False)  # run the animation
    if save:
        ani.save(title, writer='ffmpeg')
    plt.show()
