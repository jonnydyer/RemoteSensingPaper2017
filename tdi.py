#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:28:10 2017

@author: jonny
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc

rc('figure', figsize=(3.5,3))
rc('legend', fontsize='x-small')
rc('font', family='serif')

k_pe = 2e5
QE = 0.4
Q = 1.
alpha = 0.2
L_sat = 200.      # W/m^2-sr-um
lambda_c = 600e-9
delta_lambda = 0.5    # um
Vs = 7e3
Ne_well = 20000
KAI_FPS = 4
KAI_H = 6400 / 4.

def mtf_tdi_smear(phi, N, S, rho):
    return np.sinc(np.pi * rho * (1. / phi + S / N / phi)) * \
        np.sinc(np.pi * S * rho)

def tdi_equiv_linear_smear(phi, N, S, rho):
    return S + 1./phi + S / N / phi
    
def N_TDI(Q, Vs, GSD, Ne_FWC):
    LR = Vs / GSD
    return Ne_FWC * Q**2 * LR / k_pe / L_sat / delta_lambda / QE

GSD = [1., 3.]
Q = np.linspace(0.5, 1.5)
N = np.zeros((len(Q), len(GSD)))

for i, g in enumerate(GSD):
    N[:,i] = N_TDI(Q, Vs, g, Ne_well)

f = plt.figure()
f.set_tight_layout(True)
plt.plot(Q, N, linewidth=1.5)
legend = ['$N_{e^-}^{FWC}=%dke^-$, %.2f m GSD' % (Ne_well*1e-3, g) for g in GSD]

for i, g in enumerate(GSD):
    N[:,i] = N_TDI(Q, Vs, g, Ne_well*4)

plt.plot(Q, N, linewidth=1.5, linestyle='--')
legend += ['$N_{e^-}^{FWC}=%dke^-$, %.2f m GSD' % (Ne_well*4e-3, g) for g in GSD]
plt.legend(legend, fontsize='small')

plt.xlabel('Q')
plt.ylabel('$N_{TDI}$')
plt.ylim(1, 150)
plt.title('Optimal TDI stages for $V_s$ = %.1f km/s\n QE=%.2f, $L_{sat}=%d$ W/m$^2$-sr-um' %\
          (Vs/1e3, QE, L_sat), fontsize='small')
f.savefig('figures/N_tdi.pgf')

def LR_req(GSD, OS):
    return Vs * OS / GSD

GSD = np.linspace(0.5, 3.0)
OS = [1, 4, 8]

LR = np.zeros([len(GSD), len(OS)])

for i, o in enumerate(OS):
    LR[:, i] = Vs / GSD * o

f = plt.figure()
f.set_tight_layout(True)
plt.plot(GSD, LR/1e3, linewidth=1.5)
#plt.hlines(KAI_FPS*KAI_H/1e3, min(GSD), max(GSD), linestyle='--', color='r', linewidth=1.5)
plt.legend(['%dx, $N_{e-}^{FWC}$ = %d ke-' % (o, o*Ne_well/1e3) for o in OS],# + ['KAI-29050 @ 4 FPS'], 
            fontsize='small')
plt.xlabel('GSD(m)')
plt.ylabel('Linerate (klines / sec)')
plt.title('Linerate required for equivalent well depth', fontsize='small')
f.savefig('figures/LR_req.pgf')

plt.show()
