from math import *
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from optimal_Q import giqe5, rer_approx

rc('figure', figsize=(5.5,3.5))
rc('legend', fontsize='small')
rc('font', family='serif')

def alpha_eff(D, lam, Q):
    return lam / D * (1. + 1. / Q**1.35) ** (1. / 1.35)

Q = np.linspace(0.5, 2)
D = 0.35
lam = 550e-9
h = 500e3
snr0 = [15, 30] 

alpha = alpha_eff(D, lam, Q)
#GRD = h * alpha
GRD = 1. * lam / D * h
GSD = h * lam / Q / D

f, ax1 = plt.subplots()
for snr in snr0:
    GRD1 = 10**((4.4 - giqe5(h * lam / D / Q, 
                          rer_approx(Q), snr/Q, 1.))/3.32)

    ax1.plot(Q, GRD1, label=r'$SNR_{\Delta \rho}^{Q=1}=%d$' % snr)

titl = "Effective resolution," + \
r"$GRD_{eff}$" + " according to the GIQE\n model for D = %.2f m, \
alt=%d km, " % (D, h/1e3) + \
r"$\lambda$ = %d nm" % (lam * 1e9)
ax1.plot([0.5, 2], [GRD, GRD], label=r'$\frac{\lambda R}{D}$')
#ax1.plot(Q, GSD, label=r'$GSD$')
ax1.set_xlabel(r'$\frac{\lambda F^\#}{p_{px}}$')
ax1.set_ylabel('(m)')
ax1.set_ylim(0.35, 2.)
ax1.set_xlim(0.5, 2.0)

ax1.plot(Q, GSD, label='GSD')
plt.title(titl, fontsize='medium')
plt.legend()
f.set_tight_layout(True)
f.savefig('figures/resolution_q.pgf')
plt.show()
