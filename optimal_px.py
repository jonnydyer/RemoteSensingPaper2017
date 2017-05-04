#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 07:42:26 2017

@author: jonny
"""
from math import *
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

rc('figure', figsize=(3.5,2.5))
rc('legend', fontsize='x-small')

alpha15 = 0.15
alpha8 = 0.08
lambda_c = 600e-9
delta_lambda = 0.5    # um
SNR_gsd_typ = 50      # SNR / m
gsd0 = 1.
Fnum0 = 10.
Q0 = 1.0
p0 = Fnum0 * lambda_c / Q0
p = np.linspace(2, 8) * 1e-6

Ne_well0 = (SNR_gsd_typ / (sqrt(alpha15) - alpha8/sqrt(alpha15)))**2
print "Pixel well req'd: %.2f ke / um^2" % (Ne_well0 / (p0*1e6)**2 / 1e3)

Fnum = np.linspace(2, 12)
Q_opt = np.array([0.8, 1.0, 1.5], ndmin=2)

p_opt = Fnum * lambda_c / Q_opt.T
Ne_well = Ne_well0 * (p / p0)**2

f = plt.figure(1)
plt.plot(Fnum, p_opt.T*1e6)
plt.legend([r'Q = %.1f' % q for q in Q_opt.T])
plt.xlabel(r'$F^{\#}$')
plt.ylabel('Pixel Pitch (um)')
plt.title(r'$\lambda=%d$ nm' % (lambda_c * 1e9))
f.set_tight_layout(True)
f.savefig('figures/Q_pix_size.pgf')

f = plt.figure(2)
plt.plot(p*1e6, Ne_well/1e3, label='1x Oversample')
plt.plot(p*1e6, Ne_well/1e3/sqrt(2), label='2x Oversample')
plt.plot(p*1e6, Ne_well/1e3/sqrt(8), label='8x Oversample')
plt.legend()
plt.ylabel(r'$N_{e^-}^{well}$ (ke-)')
plt.xlabel('Pixel Pitch (um)')
titl = r'$\lambda=%.0f nm, SNR_{\Delta \rho}^{15 -8}/m = %.1f$, ' % (lambda_c * 1e9, 
                                                        SNR_gsd_typ)
titl += r'$\alpha=%.2f$' % alpha
print titl
plt.title(titl,fontsize='small')
f.set_tight_layout(True)
#f.savefig('figures/OS_well_depth.pgf')