#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 22:26:58 2017

@author: jonny
"""
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

rc('figure', figsize=(5.5,3.5))
rc('legend', fontsize='small')
rc('font', family='serif')

def giqe5(gsd, rer, snr, blur):
    c = [4.4, -3.32, 3.32, -5.308, -0.402, -2.92, -0.069]
    return c[0] + c[1]*np.log10(gsd) + \
        c[2]*(1. - np.exp(c[3] / snr)) * np.log10(rer) + \
            c[4]*np.log10(rer)**4 + c[5]/snr + c[6]*blur

def rer_approx(Q):
    '''
    Auelemann approximation for RER as a function of Q
    '''
    return 0.7 / (Q * (1. + 1./(Q**1.35))**(1./1.35))

def SNR_GIQE(Ne_well=10e3, Ne_rd=15.):
    return 0.07 * Ne_well / (np.sqrt(0.15 * Ne_well) + Ne_rd)

if __name__ == '__main__':

    snr0 = 30.
    rer0 = rer_approx(1)
    Fnum0 = 10.
    px_pitch0 = 6.5e-6
    d_ap0 = 0.35
    Ne_well0 = 30e3
    alt = 500e3
    lamda = 650e-9
    snr_exp = 1.         # 1 for stabilized systems or variable blur, 1.5 for constant blur
    gsd0 = px_pitch0 / (Fnum0*0.35) * alt
    NIIRS0 = giqe5(gsd0, rer0, snr0, 1.0)
    Q0 = Fnum0 * lamda / px_pitch0
    
    lamda = 550e-9
    Fnum = Fnum0
    #d_ap = [0.2, 0.35, 0.5]
    Ne_well = [20e3, 28e3, 39e3]
    
    px_pitch = np.linspace(1e-6, 20e-6, 100)
    
    iq = np.zeros([len(px_pitch), len(Ne_well)])
    Q = np.zeros([len(px_pitch), len(Ne_well)])
    gsd = np.zeros([len(px_pitch), len(Ne_well)])
    
    for i, N_e in enumerate(Ne_well):
        Q[:,i] = Fnum * lamda / px_pitch
        rer = rer_approx(Q[:,i])
        f = Fnum * d_ap0
        gsd[:,i] = px_pitch / f * alt
        #rer = rer0 * gsd[:,i] / gsd0
        #snr = snr0 / (Q[:,i]/Q0)**snr_exp
        snr = SNR_GIQE(px_pitch**2 / px_pitch0**2 * N_e)             
        iq[:,i] = giqe5(gsd[:,i], rer, snr, 1.0) - NIIRS0
    
    #plt.figure(1)
    #plt.plot(px_pitch*1e6, Q)
    #plt.xlabel('Pixel Pitch (um)')
    #plt.ylabel('Q')
    #plt.legend(['Ne_well = %.2f m' % (N/1e3) for N in Ne_well])
    #
    #plt.figure(2)
    #plt.plot(px_pitch*1e6, gsd)
    #plt.xlabel('Pixel Pitch (um)')
    #plt.ylabel('GSD (m)')
    #plt.legend(['Ne_well = %.2f ke-' % (N/1e3) for N in Ne_well])
    #
    #plt.figure(3)
    #plt.plot(Q, gsd)
    #plt.xlabel('Q')
    #plt.ylabel('GSD (m)')
    #plt.legend(['Ne_well = %.2f ke-' % (N/1e3) for N in Ne_well])
    #plt.xlim(0.5, 3)
    #
    #plt.figure(4)
    #plt.plot(Q, snr)
    #plt.xlabel('Q')
    #plt.ylabel('SNR')
    #plt.legend(['Ne_well = %.2f ke-' % (N/1e3) for N in Ne_well])
    #plt.xlim(0.5, 3)
    
    f = plt.figure(5)
    plt.plot(Q, iq)
    #plt.vlines(Fnum0*lamda/px_pitch0, -0.8, 0.8, linestyle='--', 
    #           color='r', linewidth=1)
    plt.xlabel(r'Q $\left(\frac{\lambda F^\#}{p_{px}}\right)$')
    plt.ylabel('$\Delta$ NIIRS')
    plt.legend([r'$SNR_{\Delta \rho}$ = %.0f at $Q=1$' % (SNR_GIQE(N)) for N in Ne_well] +
            ['SkySat-C'], loc='lower right')
    plt.xlim(0.5, 2)
    plt.ylim(-0.5, 0.5)
    f.set_tight_layout(True)
    #plt.grid(True)
    f.savefig('figures/Q_iq.pgf')
    
    #f = plt.figure(6)
    #plt.plot(np.linspace(0.5, 4), rer_approx(np.linspace(0.5, 4)))
    #plt.xlabel('Q')
    #plt.ylabel('RER')
    
    
    plt.show()
