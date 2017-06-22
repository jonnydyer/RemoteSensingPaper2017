# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:24:34 2016

@author: jonny
"""

from math import *
import numpy as np
import poppy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp2d, RectBivariateSpline
from scipy.signal import fftconvolve

lam = 633e-9

class gs_basis_factory:
    def __init__(self, aperture):
        self.ap=aperture
    def __call__(self, outside=None, **kwargs):
        return poppy.zernike.arbitrary_basis(self.ap.sample(), **kwargs)

class telescope(object):
    def __init__(self, d_ap, zernikes, det_os=1, obsc=0.35, 
                 npix=512, ifov=1.805e-6, lam=600e-9, fov=1e-5):
        self.d_ap = d_ap
        self.det_os = det_os
        self.lam = lam
        self.osys = poppy.OpticalSystem(oversample=det_os, npix=npix)
        primary = poppy.CircularAperture(radius=d_ap/2.)
        sec = poppy.SecondaryObscuration(secondary_radius=d_ap/2*obsc, 
                                         n_supports=0)   # secondary with 

        self.wfe = poppy.ParameterizedWFE(coefficients=np.array(zernikes)*self.lam,
                                     radius=d_ap/2,
                                     basis_factory=poppy.zernike.zernike_basis) # gs_basis = gs_basis_factory(ap)

        self.ap = poppy.CompoundAnalyticOptic(opticslist=[primary,sec, self.wfe], name='Aperture')
        self.osys.add_pupil(self.ap)
        self.osys.add_detector(pixelscale=ifov*180/pi*3600, 
                               fov_arcsec=fov*180/pi*3600, 
                               oversample=det_os)  # image plane coordinates in arcseconds
    
    def psf(self):
        return self.osys.calc_psf(self.lam)[0].data
    
    def psf_cross_section(self):
        psf = self.psf()
        psf_interp = RectBivariateSpline(np.arange(psf.shape[0]), 
                                         np.arange(psf.shape[1]), psf)
        psf_i = psf_interp.ev(np.arange(psf.shape[0]/2-self.det_os*3,
                                        psf.shape[0]/2+self.det_os*3),
                              np.ones(self.det_os*6) * psf.shape[1]/2.)
        
        return psf_i
    
    def plot_pupil(self):
        self.ap.display(what='both', opd_vmax=self.lam)
    
    def get_opd(self, npix=512):
        amplitude, pixelscale = self.ap.sample(wavelength=self.lam, npix=npix,
                                            what='amplitude', 
                                            return_scale=True)
        opd, pixelscale = self.ap.sample(wavelength=self.lam, what='opd', npix=npix,
                                      return_scale=True)
        opd[np.where(amplitude == 0)] = np.nan
        return opd
    
    def get_mtf_2d(self):
        mtf = np.abs(np.fft.fftshift(np.fft.fft2(self.psf())))
        r = np.linspace(0, mtf.shape[0]/self.det_os, mtf.shape[0])
        mtf_interp = RectBivariateSpline(np.arange(mtf.shape[0]), 
                                         np.arange(mtf.shape[1]), mtf)
        return mtf_interp, mtf, r


    def get_mtf(self):
        '''
        Return MTF curve in two directions
        '''
        mtf_interp, mtf, r = self.get_mtf_2d()
        mtf_x = mtf_interp.ev(mtf.shape[0]/2 + r,
                               np.ones(mtf.shape[1])*mtf.shape[1]/2) / np.max(mtf)
        mtf_xy = mtf_interp.ev(mtf.shape[0]/2 + r/sqrt(2),
                               mtf.shape[1]/2 + r/sqrt(2)) / np.max(mtf)
        return mtf_x, mtf_xy
    
    def get_edge_response(self):
        '''
        Compute and return edge response
        '''
        psf = self.psf()
        edge = np.zeros([psf.shape[0]*10, psf.shape[1]*10])
        edge[:, psf.shape[1]*5:] = 1
        print psf.shape, edge.shape
        er = fftconvolve(edge, psf)
        er -= np.min(er)
        er /= np.max(er)
        return er
    
    def get_edge_response1d(self):
        er = self.get_edge_response()
        er1d = er[er.shape[0]/2, :]
        xdim = np.linspace(-len(er1d)/ 2. / self.det_os, 
                           len(er1d)/ 2. / self.det_os, len(er1d))
        return xdim, er1d
    
    def get_rer(self):
        xdim, er1d = self.get_edge_response1d()
        rer = np.interp(0.5, xdim, er1d) - np.interp(-0.5, xdim, er1d)
        return rer

if __name__ == "__main__":

    ts_nom = telescope(0.35, 
                        [0.0, # piston
                        0.0, # tip
                        0.0, # tilt
                        0.15, # focus
                        0.03, # astig x
                        0.0, # astig y
                        0.00, # coma x
                        0.0, # coma y
                        0.0, # astig 2 x
                        0.0, # astig 2 y
                        -0.05],  # spherical
                      det_os=9)
    
    ts_diff = telescope(0.35, 
                        [0.0, # piston
                        0.0, # tip
                        0.0, # tilt
                        0.0, # focus
                        0.0, # astig x
                        0.0, # astig y
                        0.0, # coma x
                        0.0, # coma y
                        0.0, # astig 2 x
                        0.0, # astig 2 y
                        0.0],  # spherical
                      det_os=9)
    
    cmap = 'jet'
    opd = ts_nom.get_opd(npix=512)
    f = plt.figure(figsize=(9,5))
    
    ax = plt.subplot(221, projection='3d')
    plt.title('Wavefront Error (WFE)')
    X = np.linspace(-ts_nom.d_ap/2, ts_nom.d_ap/2, opd.shape[0])
    Y = np.linspace(-ts_nom.d_ap/2, ts_nom.d_ap/2, opd.shape[1])
    X, Y = np.meshgrid(X, Y)
    #surf = ax.plot_surface(X,Y, opd, cmap=cmap, linewidth=0.,
    #                       cstride=5, rstride=5, vmin=np.nanmin(opd), 
    #                       vmax=np.nanmax(opd), antialiased=True, alpha=1.)
    surf = ax.plot_surface(X,Y, opd, cmap=cmap, edgecolor='k', linewidth=0.5,
                           vmin=np.nanmin(opd), 
                           vmax=np.nanmax(opd), antialiased=True,
                           cstride=8, rstride=8, alpha=.5)
    ax.set_zlim(-lam/2, lam/2)
    ax.set_frame_on(False)
    ax.set_axis_off()
    ax.autoscale_view(tight=True)
    ax.view_init(elev=30., azim=0.)
    ax.dist = 5.2
    
    plt.subplot(223)
    plt.imshow(opd, cmap=cmap)
    plt.grid(False)
    plt.axis('off')
                    
    plt.subplot(222)
    plt.title('Point Spread Function (PSF)')  
    plt.imshow(ts_nom.psf(), cmap=cmap)#, vmin=0, vmax=0.1)
    plt.grid(False)
    plt.axis('off')
    
    plt.subplot(224)
    plt.title('Modulation Transfer Function (MTF)')  
    mtf_x, mtf_xy = ts_nom.get_mtf()
    plt.plot(np.linspace(0, 1, len(mtf_x)), mtf_x, label='X')
    plt.plot(np.linspace(0, 1, len(mtf_xy)), mtf_xy, label='45deg')
    plt.ylim(0, 1)
    plt.xlim(0,1)
    plt.legend(loc='best', fontsize='x-small')
    
    #f.savefig('out.png')
    
    f = plt.figure()
    rms_wfe = np.sqrt(np.mean(opd[~np.isnan(opd)]**2))/lam
    print rms_wfe
    psf_diff = ts_diff.psf_cross_section()
    psf_nom = ts_nom.psf_cross_section()
    plt.plot(psf_diff / np.max(psf_diff), label='Diffraction only')
    plt.plot(psf_nom/np.max(psf_diff), label='Aberrated')
    plt.legend(loc='best', fontsize='small')
    plt.ylim(0,1)
    plt.title(r'WFE = %.3f $\lambda$ RMS' % rms_wfe)
    
    #f.savefig('strehl.png')
    
#    Q = np.linspace(0.5, 2, 10)
#    rer = np.zeros_like(Q)
#    for i, qi in enumerate(Q):
#        ts = telescope(0.35, 
#                        [0.0, # piston
#                        0.0, # tip
#                        0.0, # tilt
#                        0.0, # focus
#                        0.0, # astig x
#                        0.0, # astig y
#                        0.0, # coma x
#                        0.0, # coma y
#                        0.0, # astig 2 x
#                        0.0, # astig 2 y
#                        0.0],  # spherical
#                         ifov = lam / 0.35 / qi,
#                         lam = lam,
#                         obsc=0.0,
#                      det_os=15)
#        rer[i] = ts.get_rer()
#    rer_approx = 1.2 / (Q * (1. + 1./(Q**1.35))**(1./1.35))
#    plt.figure()
#    plt.plot(Q, rer, label='exact')
#    plt.plot(Q, rer_approx, label='approx')
#    plt.legend()
    
    plt.show()
