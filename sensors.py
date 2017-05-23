#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 09:30:14 2017

@author: jonny
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc
import pandas as pd

rc('figure', figsize=(3.5,3))
rc('legend', fontsize='small')
rc('font', family='serif')
rc('font', size=7)

sensors = pd.read_csv('figures/sensors.csv')

sensors['kpi1'] = sensors['FWC'] * sensors['Width'] * sensors['Height'] * \
                   sensors['FPS'] / \
                   sensors['Pixel Size']

ccds = sensors.loc[sensors['Type'] == 'CCD']
cmos = sensors.loc[sensors['Type'] == 'CMOS']

f = plt.figure(1)
plt.plot(np.linspace(2,10), 0.7*np.linspace(2,10)**2, 'r--',
         label=r'$N_{e^-}^{FWC} \propto p_{px}^2$')
plt.plot(ccds['Pixel Size'], ccds['FWC']/1e3, 'b^', label='CCD')
plt.plot(cmos['Pixel Size'], cmos['FWC']/1e3, 'go', label='CMOS')
plt.xlabel(r'$p_{px}$ (um)')
plt.ylabel(r'$N_{e^-}^{FWC} (ke-)$')
plt.xlim(2, 8)
plt.ylim(0, 50)
plt.tight_layout()
plt.legend()
f.savefig('figures/p_fwc.pgf')

f = plt.figure(2)
plt.semilogy(ccds['Pixel Size'], ccds['kpi1'], 'b^', label='CCD')
plt.semilogy(cmos['Pixel Size'], cmos['kpi1'], 'go', label='CMOS')
plt.xlabel(r'$p_{px}$ (um)')
plt.ylabel(r'$\psi_{px} = \frac{h_{px}N_{e^-}^{FWC} LR}{p_{px}}$')
plt.xlim(2, 8)
plt.tight_layout()
plt.legend()
f.savefig('figures/p_kpi.pgf')

plt.show()
