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

rc('figure', figsize=(5.5,3.5))
rc('legend', fontsize='small')
rc('font', family='serif')
rc('font', size=9)

kpe = 2e5
QE = 0.4
Lsat = 200
Nbl = 1.
Q = 1.
Ne_fwc = 20e3

Vg = np.linspace(0.01, 10) * 1e3
GSD = np.linspace(0.5, 10)

frontier = kpe * Lsat * QE * Nbl * np.tile(GSD, [len(Vg), 1]).T / Q**2 / Vg / Ne_fwc

f = plt.figure(1)
CS = plt.contourf(Vg/1e3, GSD, frontier, [0, 1, 1e4], 
                  colors=['white', 'grey', 'black'])
plt.ylabel('GSD (m)')
plt.xlabel(r'$V_{gnd}$ (km/s)')
plt.xlim(0, 10)
plt.ylim(0.5,10)
plt.grid(False)
plt.title(r'$\eta_{ph}$ regime for $k_{pe}=2.5 \times 10^5$, $QE=%.1f$, ' % (QE) + \
          '\n' + \
          r'$L_{sat}=%d$, $N_{bl}=%.1f$, $Q=%.1f$, $N_e^{FWC}$ = %d ke-' % (Lsat, Nbl, Q, Ne_fwc/1e3))
plt.text(6, 2, 'Blur-limited', fontdict={'color' : 'k', 'size' : 14})
plt.text(1.5, 8, 'FWC-limited', fontdict={'color' : 'w', 'size' : 14})
plt.tight_layout()
f.savefig('figures/blur_fwc_regime.pgf')
plt.show()
