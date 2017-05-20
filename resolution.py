from math import *
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from optimal_Q import giqe5, rer_approx

rc('figure', figsize=(3.5,2.5))
rc('legend', fontsize='x-small')
rc('font', family='serif')

def alpha_eff(D, lam, Q):
    return lam / D * (1. + 1. / Q**1.35) ** (1. / 1.35)

Q = np.linspace(0.5, 2)
D = 0.35
lam = 550e-9
h = 500e3
snr0 = 21.

alpha = alpha_eff(D, lam, Q)
GRD = h * alpha
GSD = h * lam / Q / D
GRD1 = 10**((4.4 - giqe5(h * lam / D / Q, 
                          rer_approx(Q), snr0/Q, 1.))/3.32)

titl = "Effective resolution according to the GIQE\n model for D = %.2f m, \
alt=%d km, " % (D, h/1e3) + \
r"$\lambda$ = %d nm" % (lam * 1e9)

f, ax1 = plt.subplots()
ax1.plot(Q, GRD, label='GRD (approx)')
ax1.plot(Q, GRD1, label='GRD (GIQE-5)')
ax1.set_xlabel(r'Q$\left[\frac{\lambda F^\#}{p_{px}}\right]$')
ax1.set_ylabel('(m)')
ax1.set_ylim(0.35, 2)
ax1.set_xlim(0.5, 2.0)

ax1.plot(Q, GSD, label='GSD')
plt.title(titl, fontsize='small')
plt.legend(loc='best')
plt.tight_layout()
plt.gcf().savefig('figures/resolution_q.pgf')
plt.show()
