# -*- coding: utf-8 -*-
"""
Code to make a plot of SNR vs reflectance or something like that
Currently uses a fake image histogram; would look better (maybe?) with real
data. I'll hunt around a bit...
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('figure', figsize=(5.5,3.5))
#rc('figure', figsize=(9,7))
rc('legend', fontsize='xx-small')
rc('font', family='serif')
rc('font', size=9)

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
        #print(i,t)
        this_nelec = lum * t
        unsat = np.where(this_nelec<e_fwc)
        nelec[unsat] += this_nelec[unsat]
    snr = nelec / (np.sqrt(nelec + e_sigma * np.sqrt(len(exptime_arr))))
    
    return(snr, lum)


def letter_symbol(axis,x,y,letter,color="green",markersize=13,lettercolor="white"):
    xoff = 0.001*markersize*np.diff(axis.get_xlim())[0] # Klugy empirical shift
    yoff = 0.0013*markersize*np.diff(axis.get_ylim())[0]
    font = {'family': 'sans-serif',
        'color':  lettercolor,
        'weight': 'bold',
        'size': 8,
        }
    if hasattr(x, "__len__"):    # Ugly to handle arrays and scalar values
        for i in np.arange(len(x)):
            axis.plot(x[i],y[i],marker="o",color=color,markersize=markersize)
            axis.text(x[i]-xoff,y[i]-yoff,letter[i], fontdict=font)
    else:            
        axis.plot(x,y,marker="o",color=color,markersize=markersize)
        axis.text(x-xoff,y-yoff,letter,fontdict=font)
    return


def plot_hdr():
    hist, xx = fakehist()
    snr_tdi, lum = exptimes2snr(6,e_fwc = 50e3)
    snr_hdr, lum = exptimes2snr(np.array([1,2,4,8])*0.9, e_fwc=15e3)

    fig1,ax = plt.subplots()
    ax.plot(xx, snr_tdi, label = "TDI E_FWC=50k")
    ax.plot(xx, snr_hdr, label = "HDR E_FWC=20k")
    ax.legend(loc=9)
    ax.set_xlabel("Top of Atmosphere Reflectance")
    ax.set_ylabel("SNR")
    ax.set_title("SNR in different modes")
    letter_symbol(ax,[15,68,90],[110,160,170],["A","B","C"])
    ax2 = ax.twinx()
    ax2.plot(xx, hist/np.max(hist), 'y')
    ax2.fill_between(xx,np.zeros_like(hist), hist/np.max(hist), color="y", alpha = 0.2)
    ax2.set_ylabel("Histogram of pixel values (normalized to peak)", color="y", 
                   fontsize=7)
    ax2.tick_params("y", colors="y")
    fig1.set_tight_layout(True)

    # This plot would look better if we weren't using a fake histogram --
    # we'd be able to make a real histogram that isn't subject to funny
    # binning artifacts
    fig2,ax = plt.subplots()
    snrmax = np.ceil(np.max(np.concatenate((snr_tdi,snr_hdr))))
    snrnbin = 20
    snrstep = snrmax / snrnbin
    snr_x = np.linspace(0,snrmax+snrstep,snrnbin+1)
    tdi_snr_hist = np.zeros_like(snr_x)
    hdr_snr_hist = np.zeros_like(snr_x)
    for i in xx:
        tdi_snr_hist[np.int(np.round(snr_tdi[i]/snrstep-0.5))] += hist[i]
        hdr_snr_hist[np.int(np.round(snr_hdr[i]/snrstep-0.5))] += hist[i]
    snrnorm = np.ceil(np.max(np.concatenate((tdi_snr_hist, hdr_snr_hist))))
    ax.step(snr_x, tdi_snr_hist/snrnorm, label = "TDI", where="mid")
    ax.step(snr_x, hdr_snr_hist/snrnorm, label = "HDR", where="mid")
    ax.legend()
    ax.set_xlabel("SNR")
    ax.set_ylabel("Normalized Histogram of pixel SNR")
    letter_symbol(ax, [0,140,205], [6,15,3]/snrnorm, ["C", "A", "B"])
    #print(np.sum(tdi_snr_hist))
    #print(np.sum(hdr_snr_hist))
    fig2.set_tight_layout(True)
    
    return(fig1,fig2)    

if __name__ == "__main__":
    fig1,fig2 = plot_hdr()
    fig1.savefig("figures/snr_vs_L.pdf")
    fig2.savefig("figures/hist_vs_snr.pdf")
    plt.show()
