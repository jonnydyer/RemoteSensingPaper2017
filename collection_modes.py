#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 19:20:26 2017

@author: jonny
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc

rc('figure', figsize=(5.5,3.0))
rc('legend', fontsize='small')
rc('font', family='serif')

fc = (0.7,0.7,1.0)
ec = (0.2, 0.2, 0.2)
hc = (0.9, 0.6, 0.6)
lw = 0.5
alpha = 1.0

def config_axis(ax):
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def step_stare_patches(ax, ti, tf, N_it, vg, gsd, N_patch=5):
    y0 = (N_it-1)
    N_oversamp = 0
    for i in xrange(N_patch):
        for j in xrange(N_it):
            if abs(y0 - i * tf * vg - j * gsd) <= 0.5*gsd:
                c = hc
                N_oversamp += 1
            else:
                c = fc
            p = patches.Polygon([
                    [i * tf, i * tf * vg + j * gsd],
                    [i * tf + ti, (i * tf + ti) * vg + j * gsd],
                    [i * tf + ti, (i * tf + ti) * vg + (j+1) * gsd],
                    [i * tf, i * tf * vg + (j+1) * gsd],
                    ], closed=True, fc=c, ec=ec, alpha=alpha, linewidth=lw)
            ax.add_patch(p)
#    ax.plot([0,(N_patch-1) * tf + ti], [0, vg * ((N_patch-1) * tf + ti)],
#             'g--', linewidth=lw)
#    ax.plot([0,(N_patch-1) * tf + ti],
#            [N_it * gsd, vg * ((N_patch-1) * tf + ti) + N_it * gsd],
#            'g--', linewidth=lw)
    ax.set_xlim(0, N_patch * tf)
    ax.set_ylim(0, N_it * gsd + N_patch * tf * vg)
    titl = 'Step Stare $t_i(eff) = %.1f ms$\n$t_f = %d ms$, ' % \
                            (ti*N_oversamp*1e3, tf*1e3)
    titl += '$v_g = %.1f km/s$, $GSD = %.1f m$, ' % (vg/1e3, gsd)
    titl += '$t_i = %.1f ms$' % (ti*1e3)
    ax.set_title(titl, fontsize='small')
    ax.set_xticks([])
    config_axis(ax)

def stab_step_stare_patches(ax, ti, tf, N_it, vg, gsd, N_patch=5):
    y0 = (N_it-1) * gsd
    N_oversamp = 0
    for i in xrange(N_patch):
        for j in xrange(N_it):
            if abs(y0 - i * tf * vg - j * gsd) <= 0.5*gsd:
                c = hc
                N_oversamp += 1
            else:
                c = fc
            p = patches.Polygon([
                    [i * tf, i * tf * vg + j * gsd],
                    [i * tf + ti, i * tf * vg + j * gsd],
                    [i * tf + ti, i * tf * vg + (j+1) * gsd],
                    [i * tf, i * tf * vg + (j+1) * gsd],
                    ], closed=True, fc=c, ec=ec, alpha=alpha, linewidth=lw)
            ax.add_patch(p)
#    ax.plot([0,(N_patch-1) * tf + ti], [0, vg * ((N_patch-1) * tf + ti)],
#             'g--', linewidth=lw)
#    ax.plot([0,(N_patch-1) * tf + ti],
#            [N_it * gsd, vg * ((N_patch-1) * tf + ti) + N_it * gsd],
#            'g--', linewidth=lw)
    ax.set_xlim(0, N_patch * tf)
    ax.set_ylim(0, N_it * gsd + N_patch * tf * vg)

    titl = 'Stab. Step Stare $t_i(eff) = %.1f ms$\n$t_f = %d ms$, ' % \
                                   (ti*N_oversamp*1e3, tf*1e3)
    titl += '$v_g = %.1f km/s$, $GSD = %.1f m$, ' % (vg/1e3, gsd)
    titl += '$t_i = %.1f ms$' % (ti*1e3)
    ax.set_title(titl, fontsize='small')
    ax.set_xticks([])
    config_axis(ax)

def tdi_patches(ax, N_tdi, k_cc, vg, gsd, N_patch=5, color_patch=3):
    Ndot_line = vg / gsd
    ti = 1. / Ndot_line
    for i in xrange(N_patch):
        for j in xrange(N_tdi):
            if i == color_patch:
                c = hc
            else:
                c = fc
            p = patches.Polygon([
                    [(i+j) * ti, i * gsd],
                    [i * ti + (j+1) * ti, (i + k_cc) * gsd],
                    [i * ti + (j+1) * ti, (i + 1 + k_cc) * gsd],
                    [(i+j) * ti, (i+1) * gsd],
                    ], closed=True, fc=c, ec=ec, alpha=alpha, linewidth=lw)
            ax.add_patch(p)

    #ax.plot([0,N_patch * ti], [0, vg * N_patch * ti], 'g--', linewidth=lw)
#    ax.plot([N_tdi * ti,(N_patch + N_tdi) * ti],
#            [0, vg * N_patch * ti], 'g--', linewidth=lw)
    ax.set_xlim(0, (N_patch + N_tdi) * ti)
    ax.set_ylim(0, vg * N_patch * ti)
    titl = 'TDI $t_i (eff) = %.1f ms$\n$N_{tdi} = %d$, ' % (ti*N_tdi*1e3, N_tdi)
    titl += '$v_g = %.1f km/s$, $GSD = %.1f m$' % (vg/1e3, gsd)
    ax.set_title(titl, fontsize='small')
    ax.set_xticks([])
    config_axis(ax)

FPS = 300
f = plt.figure()
ax1 = plt.subplot()
step_stare_patches(ax1, ti=1./3.5e3, tf=1./FPS, N_it=36, vg=3.5e3, gsd=1, N_patch=4)
ax1.set_xlabel('Time')
ax1.set_ylabel('In-Track Dimension')
f.set_tight_layout(True)
f.savefig('figures/step_stare.pgf')

f = plt.figure()
ax1 = plt.subplot()
step_stare_patches(ax1, ti=1./3.5e3, tf=1./30, N_it=500, vg=3.5e3, gsd=1, N_patch=8)
ax1.set_xlabel('Time')
ax1.set_ylabel('In-Track Dimension')
f.set_tight_layout(True)
f.savefig('figures/step_stare_real.pgf')

f = plt.figure()
ax2 = plt.subplot()
stab_step_stare_patches(ax2, ti=1.4e-3, tf=1./FPS, N_it=36, vg=3.5e3, gsd=1, N_patch=4)
ax2.set_xlabel('Time')
ax2.set_ylabel('In-Track Dimension')
f.set_tight_layout(True)
f.savefig('figures/stab_step_stare.pgf')

f = plt.figure()
ax3 = plt.subplot()
N_tdi = 12
tdi_patches(ax3, N_tdi=N_tdi, k_cc=2.0/N_tdi, vg=3.5e3, gsd=1.,
            N_patch=32, color_patch=17)
ax3.set_xlabel('Time')
ax3.set_ylabel('In-Track Dimension')
f.set_tight_layout(True)
f.savefig('figures/tdi.pgf')

plt.show()
