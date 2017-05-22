#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:43:55 2017

@author: jonny
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc

rc('figure', figsize=(3.5,3))
rc('legend', fontsize='small')
rc('font', family='serif')
rc('font', size=7)

kpe = 2e5
QE = 0.4
Lsat = 200
Nbl = 1.
Q = 1.
Ne_fwc = 25e3

Vg = np.linspace(1, 10) * 1e3
GSD = np.linspace(0.5, 10)

frontier = kpe * Lsat * QE * Nbl * np.tile(GSD, [len(Vg), 1]).T / Q**2 / Vg / Ne_fwc

f = plt.figure(1)
plt.contourf(GSD, Vg/1e3, frontier, [0, 1, 10], colors=['white', 'grey', 'black'])
plt.xlabel('GSD (m)')
plt.ylabel(r'$V_{gnd}$ (km/s)')
plt.xlim(0.5, 10)
plt.grid(False)
plt.title(r'$\eta_{ph}$ regime for $k_{pe}=2.5 \times 10^5$, $QE=%.1f$, ' % (QE) + \
          '\n' + \
          r'$L_{sat}=%d$, $N_{bl}=%.1f$, $Q=%.1f$, $N_e^{FWC}$ = %d ke-' % (Lsat, Nbl, Q, Ne_fwc/1e3))
plt.text(5, 5, 'FWC-limited', fontdict={'color' : 'k', 'size' : 12})
plt.text(1, 8, 'Blur-limited', fontdict={'color' : 'w', 'size' : 12})
plt.tight_layout()
f.savefig('figures/blur_fwc_regime.pgf')
plt.show()