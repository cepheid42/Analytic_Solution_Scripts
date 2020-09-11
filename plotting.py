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
