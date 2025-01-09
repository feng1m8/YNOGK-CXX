import numpy as np
import matplotlib.pyplot as plt
import ynogk


def thindisk(a_spin, robs, incl, rout, m):
    pt = ynogk.ptcl(a_spin, robs, np.cos(np.deg2rad(incl)), np.sin(np.deg2rad(incl)))
    pt.metric()
    pt.center_of_image()

    alpha = np.linspace(-40, 40, m) + pt.alphac
    beta = np.linspace(-15, 15, m) + pt.betac
    alpha, beta = np.meshgrid(alpha, beta)
    alpha, beta = alpha.flatten(), beta.flatten()

    gobs = np.full(m**2, np.nan)

    rms = ynogk.rms(pt.a_spin)
    for i in range(m**2):
        pt.lambdaq(alpha[i], beta[i])

        if pt.pemdisk_all(0, rout, rms) >= 0:
            omegak = 1 / (pt.a_spin + pt.r_p**1.5)
            omegas, expnu, exppsi, _, _ = ynogk.metric(pt.r_p, 0, 1, pt.a_spin)
            vphi = exppsi / expnu * (omegak - omegas)
            gamma = 1 / np.sqrt(1 - vphi**2)
            gobs[i] = (1 - pt.omega * pt.lamda) / pt.expnu / pt.f1234[0] / (1 - omegak * pt.lamda) / (gamma / expnu)

    return -alpha.reshape(m, m), -beta.reshape(m, m), gobs.reshape(m, m)


def plot(save):
    alpha, beta, gobs = thindisk(0.998, 40, 86, 22, 1200)

    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$\beta$')
    plt.xlim(-30, 30)
    plt.ylim(-13, 15)
    plt.gca().set_aspect('equal')

    plt.pcolormesh(alpha, beta, gobs, cmap='Purples_r', rasterized=True)
    plt.colorbar(location='top', label='$g$')

    if save:
        plt.savefig(f'figure/thindisk.pdf')
        plt.savefig(f'figure/thindisk.png')
    else:
        plt.show()


if __name__ == '__main__':
    # plot(True)
    plot(False)
