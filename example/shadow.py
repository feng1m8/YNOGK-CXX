import numpy as np
import matplotlib.pyplot as plt
import ynogk


def shadow(a_spin, robs, m):
    pt = ynogk.ptcl(a_spin, robs, 0, 1)
    pt.center_of_image()

    alpha = np.linspace(-8, 4, m) + pt.alphac
    beta = np.linspace(-6, 6, m) + pt.betac
    alpha, beta = np.meshgrid(alpha, beta)
    alpha, beta = alpha.flatten(), beta.flatten()

    sigma = np.full(m**2, np.nan)

    for i in range(m**2):
        pt.lambdaq(alpha[i], beta[i])
        pt.rtp()

        if pt.r_tp1 <= pt.rhorizon:
            p_end = pt.r2p(pt.rhorizon, 0, 0)
            pt.ynogk(p_end)
            sigma[i] = pt.sigma_p
        else:
            p_end = pt.r2p(pt.robs, 1, 0)
            pt.ynogk(p_end)
            sigma[i] = 0.5 * pt.sigma_p

    return -alpha.reshape(m, m), -beta.reshape(m, m), sigma.reshape(m, m)


def plot(save):
    alpha, beta, sigma = shadow(0.998, 1e6, 800)

    plt.gcf().set_figheight(5)

    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$\beta$')
    plt.xlim(-4, 8)
    plt.ylim(-6, 6)
    plt.gca().set_aspect('equal')

    plt.pcolormesh(alpha, beta, sigma, cmap='Blues', rasterized=True)
    plt.colorbar(location='top', label=r'$\sigma$').ax.ticklabel_format(useMathText=True)

    if save:
        plt.savefig(f'figure/shadow.pdf')
        plt.savefig(f'figure/shadow.png')
    else:
        plt.show()


if __name__ == '__main__':
    # plot(True)
    plot(False)
