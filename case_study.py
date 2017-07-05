#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:27:22 2017

@author: jonny
"""
from math import *
import numpy as np
from scipy.optimize import minimize_scalar, brentq
from matplotlib import pyplot as plt
from optimal_Q import giqe5, rer_approx, SNR_GIQE
from optical_aberrations import telescope
from matplotlib import rc

rc('figure', figsize=(3.5,2.5))
rc('legend', fontsize='x-small')
rc('font', family='serif')

mu_earth = 3.986004418e14     # m^3 / s^2
Re = 6385e3             # m
alt = 500e3             # m
lam = [450e-9, 800e-9]  # m
lam_c = 525e-9          # m
Lsat = 200.             # W / m^2-sr
Ltyp = 10.              # W / m^2-sr

alpha = Ltyp / Lsat
a = Re + alt
t_orb = 2 * pi * sqrt(a**3 / mu_earth)
Vgnd = 2 * pi * Re / t_orb

def find_min_d(niirs_req, Q, snr_giqe, alt, obsc=0.3):
    #def f(qi):
    def f(d):
        ts = telescope(d, 
            [0.0,   # piston
             0.0,   # tip
             0.0,   # tilt
             0.03,  # focus
             0.0,   # astig x
             0.05,  # astig y
             0.05,  # coma x
             0.0,   # coma y
             0.0,   # astig 2 x
             0.0,   # astig 2 y
             0.0], # spherical
             ifov = lam_c / d / Q,
             lam = lam_c,
             obsc=obsc,
             det_os=3)
        rer = ts.get_rer()
        #print 'For q=%.3f, rer=%.3f' % (qi, rer)
        gsd = alt * lam_c / Q / d
        niirs = giqe5(gsd, rer, snr_giqe, 1.0)
        return niirs-niirs_req
    #q_opt = minimize_scalar(f, method='bounded', bounds=[0.8, 2.0])
    dmin = brentq(f, 0.3, 0.7)
    #niirs_opt = f(dmin) + niirs_req
    gsd = alt * lam_c / Q / dmin
    return dmin, gsd

systems = {
    '50cm Resolution Point' : {
        'grd_eff' : 0.5,
        'snr_req' : 70.,
        'DR_req' : 75.,      # dB
        'Q' : 1.3,
        'x_ct' : 8e3,      # m
        'N_bands' : 5,
        'p_px' : 3.45e-6,  # m'
        'hpx_det' : 4096,
        'vpx_det' : 3000,
        'Ne_FWC_det' : 10e3,
        'rd_noise' : 5.
    },
    '1m Resolution Large Area' : {
        'grd_eff' : 1.0,
        'snr_req' : 70.,
        'DR_req' : 75.,      # dB
        'Q' : 1.0,
        'x_ct' : 18e3,       # m
        'N_bands' : 5,
        'p_px' : 3.45e-6,    # m
        'hpx_det' : 4096,
        'vpx_det' : 3000,
        'Ne_FWC_det' : 10e3,
        'rd_noise' : 5.
    }
}

for i, (n, s) in enumerate(systems.items()):
    print 'System: %s' % n
    niirs = -3.32 * np.log10(s['grd_eff']) + 4.4
    NsNe_fwc = s['snr_req']**2 / alpha
    print '\tNsNe_fwc = %.1f ke-' % (NsNe_fwc/1e3)
    snr_delta_rho = SNR_GIQE(NsNe_fwc, 5.)
    print '\tSNR (GIQE, 5e- rd.) = %.1f' % snr_delta_rho
    snr_delta_rho = SNR_GIQE(NsNe_fwc, 15.)
    print '\tSNR (GIQE, 15e- rd.) = %.1f' % snr_delta_rho
    print '\tNIIRS = %.2f' % niirs
    
    d, gsd = find_min_d(niirs, s['Q'], snr_delta_rho, alt)
    LR = Vgnd / gsd
    hpx = s['x_ct'] / gsd
    psi_px = NsNe_fwc * hpx * Vgnd / gsd * s['N_bands'] # p_px
    f = s['p_px'] / gsd * alt
    N_sensors = floor(hpx/s['hpx_det'])
    FPS_req = psi_px / s['Ne_FWC_det'] / hpx / N_sensors / s['vpx_det']
    DR_lin = 10**(s['DR_req'] / 20.)
    Ns_req = DR_lin**2 / (s['Ne_FWC_det']**2 / s['rd_noise'])
    Ns = floor(gsd / Vgnd * s['vpx_det'] * FPS_req / s['N_bands'])
    I_D_raw = Ns * np.log2(DR_lin)
    print '\tGSD: %.3f m' % (gsd)
    print '\tD_ap = %0.3f m' % d
    print '\tf = %.3f m' % f
    print '\tF# = %0.1f' % (f/d)
    print '\t%0.1fmm focal plane width' % (hpx * s['p_px'] * 1e3)
    print '\t%d Cross-track pixels required' % hpx
    print '\t%d Cross track sensors' % (N_sensors)
    print '\tRequired framerate %0.1f' % (FPS_req)
    print '\tpsi_px required: %.1e' % (psi_px / N_sensors)
    print '\tRequired Ns for DR_req : %0.1f' % (Ns_req)
    print '\tActual Ns : %d' % Ns
    print '\tI_D = %0.1f bit / pix raw' % (I_D_raw)
    
