import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import ynogk


def time(pt):
    rms = ynogk.rms(pt.a_spin)
 
    def func(alpha):
        pt.lambdaq(alpha, 0)
        pt.pemdisk(0, pt.robs, pt.rhorizon)
        return pt.r_p - rms

    sol = optimize.root_scalar(func, bracket=[6, 7])
    pt.lambdaq(sol.root, 0)
    pem = pt.pemdisk(0, pt.robs, pt.rhorizon)
    pt.ynogk(pem)
    return pt.time_p


def hotspot(a_spin, robs, incl, Rspot, m, n):
    pt = ynogk.ptcl(a_spin, robs, np.cos(np.deg2rad(incl)), np.sin(np.deg2rad(incl)))
    pt.metric()
    pt.center_of_image()

    alpha = np.linspace(-10, 10, n) + pt.alphac
    beta = np.linspace(-9, 5, n) + pt.betac
    alpha, beta = np.meshgrid(alpha, beta)
    alpha, beta = alpha.flatten(), beta.flatten()

    rms = ynogk.rms(-pt.a_spin)
    omegak = 1 / (pt.a_spin - rms**1.5)
    per = 2 * np.pi * (rms**1.5 - pt.a_spin)
    time_0 = time(pt)

    gobs = np.linspace(0.4, 1.5, m)
    tobs = np.linspace(0, 1, m + 1)
    tobs = 0.5 * (tobs[1:] + tobs[:-1])
    jobs = np.full((len(gobs), len(tobs)), 0, dtype=float)

    phis = (1.5 - 2 * tobs) * np.pi
    spot = np.sqrt(rms**2 + pt.a_spin**2) * np.array([np.cos(phis), np.sin(phis)])

    for i in range(len(alpha)):
        pt.lambdaq(alpha[i], beta[i])
        pem = pt.pemdisk(0, rms + 4 * Rspot, rms - 4 * Rspot)

        if pem >= 0:
            pt.ynogk(pem)
            emis = np.sqrt(pt.r_p**2 + pt.a_spin**2) * np.array([np.cos(pt.phi_p), np.sin(pt.phi_p)])
            omega, expnu, exppsi, _, _ = ynogk.metric(pt.r_p, 0, 1, pt.a_spin)
            vphi = exppsi / expnu * (omegak - omega)
            gamma = 1 / np.sqrt(1 - vphi**2)
            g = expnu / gamma / (1 - omegak * pt.lamda)

            dist = (spot[0] - emis[0])**2 + (spot[1] - emis[1])**2
            jx = np.exp(-0.5 * dist / Rspot**2)
            jx[dist > (4 * Rspot)**2] = 0

            gj = np.searchsorted(gobs, g)
            tj = np.searchsorted(tobs, (tobs + (pt.time_p - time_0) / per) % 1)
            jobs[gj - 1][tj - 1] += jx

    jobs = jobs / np.max(jobs)
    jobs[jobs == 0] = np.nan
    return tobs, gobs, jobs


def plot(save):
    tobs, gobs, jobs = hotspot(0, 1e8, 60, 0.5, 400, 800)

    plt.gcf().set_figheight(5)

    plt.xlabel('Observer Time [Periods]')
    plt.ylabel(r'$E_\mathrm{obs} / E_\mathrm{em}$')
    plt.xlim(0, 1)
    plt.ylim(0.4, 1.5)
    plt.gca().set_aspect(0.8)

    plt.pcolormesh(tobs, gobs, jobs, cmap='Reds_r', rasterized=True)
    plt.colorbar(location='top', label=r'$j_\nu(\mathbf{x})$ [Normalized]')

    if save:
        plt.savefig(f'figure/hotspot.pdf')
        plt.savefig(f'figure/hotspot.png')
    else:
        plt.show()


if __name__ == '__main__':
    # plot(True)
    plot(False)
