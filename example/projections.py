import numpy as np
import matplotlib.pyplot as plt
import ynogk


def projections(a_spin, robs, thetaobs, m, n):
    pt = ynogk.ptcl(a_spin, robs, np.cos(np.deg2rad(thetaobs)), np.sin(np.deg2rad(thetaobs)))
    pt.center_of_image()

    alpha0 = np.linspace(-10, 10, n) + pt.alphac
    beta0 = np.linspace(-10, 10, m) + pt.betac
    alpha, beta = np.meshgrid(alpha0, beta0)
    alpha0 = np.linspace(-10, 10, m) + pt.alphac
    beta0 = np.linspace(-10, 10, n) + pt.betac
    alpha0, beta0 = np.meshgrid(alpha0, beta0)
    alpha = np.concatenate((alpha, alpha0.T)).flatten()
    beta = np.concatenate((beta, beta0.T)).flatten()

    x = np.full(2 * m * n, np.nan)
    y = np.full(2 * m * n, np.nan)

    for i in range(2 * m * n):
        pt.lambdaq(alpha[i], beta[i])
        pem = pt.pemdisk(0, pt.robs, pt.rhorizon)

        if pem >= 0:
            pt.phi(pem)
            x[i] = np.sqrt(pt.r_p**2 + pt.a_spin**2) * np.cos(pt.phi_p)
            y[i] = np.sqrt(pt.r_p**2 + pt.a_spin**2) * np.sin(pt.phi_p)

    return x.reshape(2 * m, n), y.reshape(2 * m, n)


def plot(save):
    a_spin = 0.95
    x, y = projections(a_spin, 4e70, 60, 41, 400)

    plt.gcf().set_figheight(4)

    plt.xlabel('$-Y$ [$GM/c^2$]')
    plt.ylabel('$-X$ [$GM/c^2$]')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.gca().set_aspect('equal')

    for i in range(x.shape[0]):
        plt.plot(-y[i], -x[i])

    plt.gca().add_artist(plt.Circle((0, 0), radius=np.sqrt(2 * ynogk.rhorizon(a_spin)), zorder=2, color='white'))

    if save:
        plt.savefig(f'figure/projections.pdf')
        plt.savefig(f'figure/projections.png')
    else:
        plt.show()


if __name__ == '__main__':
    # plot(True)
    plot(False)
