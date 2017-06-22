# -*- coding: utf-8 -*-
"""
Code to make a plot of SNR vs reflectance or something like that
Currently uses a fake image histogram; would look better (maybe?) with real
data. I'll hunt around a bit...
"""
import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, mu, sig, A):
    return A * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def fakehist():
    xx = np.arange(100)

    asphalt = gaussian(xx, 6, 3, 10)
    grass = gaussian(xx,16, 5, 8)
    concrete = gaussian(xx,80,20, 1)
    
    return(asphalt + grass + concrete, xx)


def exptimes2snr(exptimes, e_fwc = 3e4, e_sigma = 15):
    """
    Given an exposure time, returns an SNR as a function of luminosity
    If exptimes is an array/list, assumes an HDR mode
    """
    exptime_arr = exptimes if hasattr(exptimes, "__len__") else np.array([exptimes])
    lum = np.arange(100)*1e2    # electrons per 0.1us of exposure
    nelec = np.zeros_like(lum)
    for i,t in enumerate(np.array(exptime_arr)):
        print(i,t)
        this_nelec = lum * t
        unsat = np.where(this_nelec<e_fwc)
        nelec[unsat] += this_nelec[unsat]
    snr = nelec / (np.sqrt(nelec + e_sigma * np.sqrt(np.array(exptime_arr)).size))
    
    return(snr, lum)

def plot_hdr():
    hist, xx = fakehist()
    snr_tdi, lum = exptimes2snr(6,e_fwc = 50e3)
    snr_hdr, lum = exptimes2snr(np.array([1,2,4,8])*0.9, e_fwc=15e3)

    fig,ax = plt.subplots(figsize=[10,7])
    ax.plot(xx, snr_tdi, label = "TDI E_FWC=50k")
    ax.plot(xx, snr_hdr, label = "HDR E_FWC=20k")
    ax.legend()
    ax.set_xlabel("Luminosity")
    ax.set_ylabel("SNR")
    ax.set_title("SNR in different modes")
    ax.plot([15,68,90], [110,160,170], 'ro',markersize=20)
    ax.text(14.5,108,"A",color="white")
    ax.text(67.5,158,"B",color="white")
    ax.text(89.5,168,"C",color="white")
    ax2 = ax.twinx()
    ax2.plot(xx, hist, 'y--')
    ax2.set_ylabel("Histogram of pixel values", color="y")
    ax2.tick_params("y", colors="y")

    # This plot would look better if we weren't using a fake histogram --
    # we'd be able to make a real histogram that isn't subject to funny
    # binning artifacts
    fig,ax = plt.subplots(figsize=[10,7])
    snrmax = np.ceil(np.max(np.concatenate((snr_tdi,snr_hdr))))
    snrnbin = 20
    snr_x = np.linspace(0,snrmax,snrnbin)
    tdi_snr_hist = np.zeros_like(snr_x)
    hdr_snr_hist = np.zeros_like(snr_x)
    for i in xx:
        tdi_snr_hist[np.int(np.round(snr_tdi[i]/snrnbin))] += hist[i]
        hdr_snr_hist[np.int(np.round(snr_hdr[i]/snrnbin))] += hist[i]
    ax.plot(snr_x, tdi_snr_hist, label = "TDI")
    ax.plot(snr_x, hdr_snr_hist, label = "HDR")
    ax.legend()
    ax.set_xlabel("SNR")
    ax.set_ylabel("Histogram of pixel values")

    print(np.sum(tdi_snr_hist))
    print(np.sum(hdr_snr_hist))
    return(fig)    