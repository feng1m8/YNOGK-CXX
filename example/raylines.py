import numpy as np
import matplotlib.pyplot as plt
import ynogk


def raylines(a_spin, m, n):
    rms = ynogk.rms(a_spin)
    omega, expnu, exppsi, _, _ = ynogk.metric(rms, 0, 1, a_spin)
    vphi = exppsi / expnu * (1 / (rms**1.5 + a_spin) - omega)

    pt = ynogk.ptcl(a_spin, rms, 0, 1, [0, 0, vphi])

    theta = np.linspace(0, np.pi, m)[1:-1]
    phi = np.linspace(0, 2 * np.pi, m - 1, endpoint=False)
    theta, phi = np.meshgrid(theta, phi)
    theta, phi = theta.flatten(), phi.flatten()
    theta = np.append(theta, [0, np.pi])
    phi = np.append(phi, [0, 0])

    pradius, ptheta = np.sin(theta) * [np.sin(phi), np.cos(phi)]
    pphi = np.cos(theta)

    x = np.full((2 + (m - 2) * (m - 1), n), np.nan)
    y = np.full((2 + (m - 2) * (m - 1), n), np.nan)

    for i in range(2 + (m - 2) * (m - 1)):
        pt.lambdaq(pradius[i], ptheta[i], pphi[i])
        pem = np.linspace(0, pt.p_total(), n)

        for j in range(n):
            pt.phi(pem[j])

            if pt.r_p > 1.1 * pt.rhorizon:
                x[i][j] = np.sqrt(pt.r_p**2 + pt.a_spin**2) * pt.sin_p * np.cos(pt.phi_p)
                y[i][j] = -np.sqrt(pt.r_p**2 + pt.a_spin**2) * pt.sin_p * np.sin(pt.phi_p)

    return x, y


def plot(save):
    a_spin = 0.9375
    x, y = raylines(a_spin, 21, 400)

    plt.gcf().set_figheight(5)

    plt.ylabel('$Y$ [$GM/c^2$]')
    plt.xlabel('$X$ [$GM/c^2$]')
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.gca().set_aspect('equal')

    for i in range(x.shape[0]):
        plt.plot(x[i], y[i])

    plt.gca().add_artist(plt.Circle((0, 0), radius=np.sqrt(2 * ynogk.rhorizon(a_spin)), zorder=2, fill=False, linewidth=3))

    if save:
        plt.savefig(f'figure/raylines.pdf')
        plt.savefig(f'figure/raylines.png')
    else:
        plt.show()


if __name__ == '__main__':
    # plot(True)
    plot(False)
