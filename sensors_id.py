#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Ripped from sensors.py
Added a couple of different scatter plots based on the same sensor database

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc
import pandas as pd

scaler=2.5
rc('figure', figsize=(3.5,3))
rc('legend', fontsize='x-small')
#rc('figure', grid=True)
rc('font', family='serif')
rc('font', size=7)
rc('lines', markersize=3.0)


def LoadSensorData():
    """
    Just read in the data and calculate some derived metrics from Jonny's
    sheet of sensor properties
    """
    sensors = pd.read_csv('figures/sensors.csv')
    # JD Note: changed this to match my updated KPI which doesn't put p_px in 
    # denominator because I don't think that made sense.
    sensors['kpi1'] = sensors['FWC'] * sensors['Width'] * sensors['Height'] * \
                       sensors['FPS']# / \
                 #      sensors['Pixel Size']
    
    bits = 12 # Could use log2(DR) but why penalize sensors with good noise characteristics?
    
    # Data intensity in bits/sec
    sensors['I_d'] = sensors['Width'] * sensors['Height'] * sensors['FPS'] * bits
    # Number of samples needed for SNR=100 at alpha=5
    sensors['n_samp'] = np.ceil(5e4/sensors['FWC']) 
    # Number of samples for each ground point for GSD=1m at h=500km (v=7000m/s)
    sensors['n_samp_gsd1'] = np.floor(sensors['Height']/(7000./1./sensors['FPS']))
    # Number of samples for each ground point for GSD=4m at h=500km (v=7000m/s)
    sensors['n_samp_gsd4'] = np.floor(sensors['Height']/(7000./4./sensors['FPS']))
    # Number of bands that can be sampled at SNR=100 with GSD=1m. Not floored/rounded...
    sensors['n_bands_gsd1'] = sensors['n_samp_gsd1']/sensors['n_samp']
    # Number of bands that can be sampled at SNR=100 with GSD=1m. Not floored/rounded...
    sensors['n_bands_gsd4'] = sensors['n_samp_gsd4']/sensors['n_samp']
    # Data intensity assuming on-board summing with shift/add
    # Could multiply by N_bands below but that's less fundamental
    sensors['I_d_sum'] = sensors['I_d']/sensors['n_samp_gsd1']
    
    return(sensors)


def SensorType(sensors):
    ccds = sensors.loc[sensors['Type'] == 'CCD']
    cmos = sensors.loc[(sensors['Type'] == 'CMOS') & (sensors['Shutter'] == 'Global')]
    cmos_rolling = sensors.loc[(sensors['Type'] == 'CMOS') & (sensors['Shutter'] == 'Rolling')]
    return(ccds, cmos, cmos_rolling)

def PlotPsiVsId(sensors):
    """
    Make some plots of 
    """
    ccds, cmos, cmos_rolling = SensorType(sensors)
    
    # Information content vs data intensity. Scatter about the linear fit is just
    # due to higher order variation in e.g. aspect ratio, bit depth, etc.
    fig1 = plt.figure(1)
    plt.loglog(ccds['I_d'], ccds['kpi1'], 'r^', label='CCD')
    plt.loglog(cmos['I_d'], cmos['kpi1'], 'go', label='CMOS (Global)')
    plt.loglog(cmos_rolling['I_d'], cmos_rolling['kpi1'], 'bs', label='CMOS (Rolling)')
    plt.xlabel(r'$BW_{det}$ [bits/sec]')
    plt.ylabel(r'$\psi_{px} = h_{px}N_{e^-}^{FWC} LR$')
    #plt.xlim(2, 8)
    #plt.ylim(0, 50)
    plt.tight_layout()
    plt.legend()
    
    # Information content vs data intensity assuming on-board summing of digital
    # oversamples. This (of course) favors CMOS detectors, shifting them to the
    # left so that they are above CCDs instead of up and to the right.
    fig2 = plt.figure(2)
    plt.loglog(ccds['I_d_sum'], ccds['kpi1'], 'r^', label='CCD')
    plt.loglog(cmos['I_d_sum'], cmos['kpi1'], 'go', label='CMOS (Global)')
    plt.loglog(cmos_rolling['I_d_sum'], cmos_rolling['kpi1'], 'bs', label='CMOS (Rolling)')
    plt.xlabel(r'$BW_{det}$ [bits/sec]')
    plt.ylabel(r'$\psi_{px} = h_{px}N_{e^-}^{FWC} LR$')
    #plt.xlim(2, 8)
    #plt.ylim(0, 50)
    plt.tight_layout()
    plt.legend()
    
    return(fig1, fig2)


def PlotNBands(sensors):
    # Klugy attempt to show how the height of a CMOS detector can be used to 
    # add more spectral bands
    font = {'family': 'sans_serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 7,
            }
    # Some interpolated grids for plotting the number of bands on a non-linear
    # Y axis
    yax_grid = np.array([0,0.5,0.6,1,1.5,1.6,2,3,4,5,6,7,10,15,20,25,30])
    sensors['yax'] = np.interp(sensors['n_bands_gsd1'], yax_grid, np.arange(len(yax_grid)))
    yax_grid4 = np.array([0,0.5,0.6,1,1.5,1.6,2,3,4,5,6,7,10,15,20,25,30,50])
    sensors['yax4'] = np.interp(sensors['n_bands_gsd4'], yax_grid4, np.arange(len(yax_grid4)))

    ccds, cmos, cmos_rolling = SensorType(sensors)

    fig1,ax = plt.subplots()
    ax.plot(ccds['Pixel Size'], ccds['yax'], 'r^', label='CCD')
    ax.plot(cmos['Pixel Size'], cmos['yax'], 'go', label='CMOS (Global)')
    ax.plot(cmos_rolling['Pixel Size'], cmos_rolling['yax'], 'bs', label='CMOS (Rolling)')
    ax.set_xlabel(r'$p_{px}$ ($\mu$m)')
    ax.set_ylabel(r'$N_{bands}$')
    tickind = [0,3,6,7,8,9,12,15]
    ax.set_yticks(tickind)
    ax.set_yticklabels(np.array(yax_grid[tickind]))
    ax.set_xlim(2, 8)
    ax.set_title("Channels available for GSD=1m, SNR=100")
    plt.tight_layout()
    ax.legend()
    labelx = 6.2
    ax.add_patch(patches.Rectangle((2,-1),6,4,facecolor='grey', alpha=0.2))
    ax.text(labelx, 0.2, "SNR<100",fontdict=font)
    ax.add_patch(patches.Rectangle((2,3),6,4,facecolor='grey', alpha=0.1))
    ax.text(labelx, 4, "Panchromatic",fontdict=font)
    #ax.add_patch(patches.Rectangle((2,6),6,3,facecolor='red', alpha=0.1))
    #ax.text(6.8, 7, "False color",fontdict=font)
    ax.add_patch(patches.Rectangle((2,7),6,5,facecolor='green', alpha=0.1))
    ax.text(labelx, 9.5, "Multispectral",fontdict=font)
    ax.add_patch(patches.Rectangle((2,12),6,5,facecolor='blue', alpha=0.1))
    ax.text(labelx, 13, "Hyperspectral",fontdict=font)
    #f.savefig('figures/p_kpi.pgf')
    
    fig2,ax = plt.subplots()
    ax.plot(ccds['Pixel Size'], ccds['yax4'], 'r^', label='CCD')
    ax.plot(cmos['Pixel Size'], cmos['yax4'], 'go', label='CMOS (Global)')
    ax.plot(cmos_rolling['Pixel Size'], cmos_rolling['yax4'], 'bs', label='CMOS (Rolling)')
    ax.set_xlabel(r'$p_{px}$ ($\mu$m)')
    ax.set_ylabel(r'$N_{bands}$')
    tickind = [0,3,6,7,8,9,12,15]
    ax.set_yticks(tickind)
    ax.set_yticklabels(np.array(yax_grid4[tickind]))
    ax.set_xlim(2, 8)
    ax.set_title("Channels available for GSD=4m, SNR=100")
    plt.tight_layout()
    ax.legend()
    labelx = 6.2
    ax.add_patch(patches.Rectangle((2,-1),6,4,facecolor='grey', alpha=0.2))
    ax.text(labelx, 0.2, "SNR<100",fontdict=font)
    ax.add_patch(patches.Rectangle((2,3),6,4,facecolor='grey', alpha=0.1))
    ax.text(labelx, 4, "Panchromatic",fontdict=font)
    #ax.add_patch(patches.Rectangle((2,6),6,3,facecolor='red', alpha=0.1))
    #ax.text(6.8, 7, "False color",fontdict=font)
    ax.add_patch(patches.Rectangle((2,7),6,5,facecolor='green', alpha=0.1))
    ax.text(labelx, 9, "Multispectral",fontdict=font)
    ax.add_patch(patches.Rectangle((2,12),6,6,facecolor='blue', alpha=0.1))
    ax.text(labelx, 15, "Hyperspectral",fontdict=font)

    return(fig1,fig2)


if __name__ == "__main__":
    sensors = LoadSensorData()
    fig1,fig2 = PlotPsiVsId(sensors)
    fig3,fig4 = PlotNBands(sensors)
    fig1.savefig("figures/psi_vs_id.pgf")
    fig2.savefig("figures/psi_vs_id_summed.pgf")
    fig3.savefig("figures/nbands_1m.pgf")
    fig4.savefig("figures/nbands_4m.pgf")
    plt.show()
