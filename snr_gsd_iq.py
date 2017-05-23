#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:59:14 2017

@author: jonny
"""

from math import *
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from optimal_Q import giqe5, rer_approx

rc('figure', figsize=(3.5,2.5))
rc('legend', fontsize='x-small')
rc('font', family='serif')

snr = np.linspace(5, 30)
gsd = np.linspace(0.2, 10)
rer0 = 0.3

niirs = giqe5(np.tile(gsd, [len(snr),1]).T, rer0, snr, 1)
grd = 10**((4.4 - niirs)/3.32)
snr2_gsd2 = snr**2 / np.tile(gsd, [len(snr),1]).T**2
                   
plt.figure()
CS = plt.contour(snr, gsd, grd)
plt.clabel(CS)

plt.figure()
CS = plt.contour(snr, gsd, snr2_gsd2)
plt.clabel(CS)