#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 07:41:26 2017

@author: jonny
"""

from math import *
import numpy as np
from matplotlib import rc, gridspec, transforms
import matplotlib.pyplot as plt
from optical_aberrations import telescope
from scipy.signal import fftconvolve
from skimage.measure import block_reduce

rc('figure', figsize=(3.5,3.5))
rc('legend', fontsize='x-small')
rc('axes', grid=False)

Q = 1.
d_ap = 0.35
lam = 600e-9
ifov = lam / d_ap / Q

SNR = (0.5, 50)

OS = 8 
img_dim = int(128*Q)
rho_max = 1.0# / Q
k = rho_max / img_dim

tel = telescope(d_ap,
        [0.0, # piston
        0.0, # tip
        0.0, # tilt
        0.1, # focus
        0.0, # astig x
        0.0, # astig y
        0.0, # coma x
        0.0, # coma y
        0.0, # astig 2 x
        0.0, # astig 2 y
        -0.0],  # spherical
      det_os=OS, ifov=ifov, obsc=0.3, fov=7.5*ifov)

rho = np.indices([img_dim*OS, img_dim*OS])[1].astype(np.float64) \
        / OS / img_dim
rho_small = np.indices([img_dim, img_dim])[1].astype(np.float64) / img_dim
#SNR_v = np.indices([img_dim, img_dim])[0].astype(np.float64) / img_dim * (SNR[1]-SNR[0]) + SNR[0]
#sigma = 1. / SNR_v
sigma = np.indices([img_dim, img_dim])[0].astype(np.float64) / img_dim / SNR[0] + 1./SNR[1]
SNR_v = 1. / sigma
mtf = np.abs(np.fft.fft2(tel.psf(), s=(int(img_dim*OS/rho_max*Q),
    int(img_dim*OS/rho_max*Q))))[:img_dim, :img_dim]

mtf_w_samp = mtf * np.sinc(rho_max * rho_small / Q).T

#img = np.random.normal(loc=mtf * np.sign(np.sin(2 * pi * rho**2 * k/4)), scale=sigma)
x = rho * img_dim / np.sqrt(Q)
img = (1+np.sin(2 * pi * x**2 * k/2))/2
#img = (np.sign(img - 0.5) + 1.0001) / 2
img = fftconvolve(img, tel.psf(), 'same')
img_sm = block_reduce(img, (OS, OS))
img = np.random.poisson(img_sm * SNR_v**2).astype(np.float64) / SNR_v**2

gs = gridspec.GridSpec(5, 5)
f = plt.figure(1)
ax = plt.subplot(gs[1:, 1:])
plt.imshow(img, cmap='gray', interpolation='nearest')
#plt.clim(0,1)
#CS = plt.contour(SNR_v * np.tile(mtf_w_samp[:,0], (img_dim,1)), [0.5, 1, 2], colors='y')
#plt.clabel(CS, inline=1, fontsize=10)
ax.set_axis_off()

rot = transforms.Affine2D().rotate_deg(270)
ax = plt.subplot(gs[1:, 0])
ax.plot(10*np.log10(SNR_v[:, 0]), transform = rot + ax.transData)
ax.axes.get_yaxis().set_visible(False)
ax.set_xticks(np.linspace(0, 20, 3))
ax.grid(True)
ax.set_title('SNR (dB)', fontsize='small')


ax = plt.subplot(gs[0, 1:])
plt.plot(mtf_w_samp[:,0])
ax.axvline(Q/2*len(mtf_w_samp[:,0]), linestyle='--', color='r')
#ax.axes.get_xaxis().set_visible(False)
ax.grid(True)
ax.set_yticks([0., 0.5, 1.0])
ax.set_title('MTF', fontsize='small')
ax.yaxis.tick_right()
ax.set_xticks([0, len(mtf_w_samp[:,0])])
ax.set_xticklabels([r'$0$', r'$\frac{%.1f}{\lambda F^\#}$' % rho_max])
ax.set_xlim(0,len(mtf_w_samp[:,0]))
ax.xaxis.tick_top()
#f.set_tight_layout(True)
#f.savefig('figures/SNR_mtf_Q1.pdf')

#plt.figure()
#plt.imshow(tel.psf(), cmap='viridis')

plt.show()
