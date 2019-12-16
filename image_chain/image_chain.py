from math import *
import numpy as np
import cairosvg
import imageio
import yaml
from colour_demosaicing import (demosaicing_CFA_Bayer_bilinear,
                                demosaicing_CFA_Bayer_Malvar2004,
                                demosaicing_CFA_Bayer_Menon2007)

# Measured in InkScape for letter E
# 20/20 vision implies resolving characters with line spacing of 1 arc-min
# The letter E has 5 such lines vertically and was measured 27.33 px high in
# Inkscape.

SNELLEN_20_20_SCALE = (5. * pi / 180. / 60.) / 27.33
ARC_MIN_IN_RAD = pi / 180. / 60.

DEMOSAIC_ALGOS = {"bilinear" : demosaicing_CFA_Bayer_bilinear,
         "malvar" : demosaicing_CFA_Bayer_Malvar2004,
         "menon" : demosaicing_CFA_Bayer_Menon2007}

def get_snellen_raster(ifov, oversample=16):
    '''
    Returns a raster snellen chart scaled such that a pixel in the image is
    1 arc-min / oversample in extent.
    '''
    cairosvg.svg2png(file_obj = open("snellen_chart.svg", "rb"),
                     write_to = "/tmp/output.png",
                     scale=SNELLEN_20_20_SCALE * oversample / ifov)
    out = imageio.imread("/tmp/output.png")
    return out.astype(np.float64) / 255.

class Detector(object):
    def __init__(self, fname="conf.yml"):
        with open(fname, 'r') as f:
            self._params = yaml.load(f)
            for k,v in self._params.items():
                # Load config params in as attributes unless
                # attribute already exists (as a property for instance)
                setattr(self, k, v)

    def apply_filter_gain_rgb(self, img):
        img[:, :, 0] *= self.detector["red_gain"]
        img[:, :, 1] *= self.detector["green_gain"]
        img[:, :, 2] *= self.detector["blue_gain"]

    def invert_filter_gain_rgb(self, img):
        img[:, :, 0] /= self.detector["red_gain"]
        img[:, :, 1] /= self.detector["green_gain"]
        img[:, :, 2] /= self.detector["blue_gain"]

    def apply_noise_rgb(self, img):
        img[:, :, 0] = np.random.poisson(img[:, :, 0] *
                                         self.detector['well_depth'])
        img[:, :, 1] = np.random.poisson(img[:, :, 1] *
                                         self.detector['well_depth'])
        img[:, :, 2] = np.random.poisson(img[:, :, 2] *
                                         self.detector['well_depth'])
        img += (np.random.random(img.shape) * self.detector['rd_noise'] -
                self.detector['rd_noise']/2)
        img /= np.max(img)

    def apply_noise_raw(self, img):
        img2 = np.random.poisson(img *
                                self.detector['well_depth']).astype(np.float64)
        img2 += (np.random.random(img2.shape) * self.detector['rd_noise'] -
                self.detector['rd_noise']/2)
        img2 /= np.max(img2)
        return img2

    def get_raw_img(self, img):

        dim_out_2 = (int(np.ceil(float(img.shape[0]) / 2.)),
                   int(np.ceil(float(img.shape[1]) / 2.)))
        img_raw = np.zeros((2 * dim_out_2[0], 2 * dim_out_2[1]))

        red = np.tile(np.array([[self.detector["red_gain"], 0.0],
                                [0.0, 0.0]]), (dim_out_2[0], dim_out_2[1]))
        green = np.tile(np.array([[0.0, self.detector["green_gain"]],
                         [self.detector["green_gain"], 0.0]]),
                        (dim_out_2[0], dim_out_2[1]))
        blue = np.tile(np.array([[0.0, 0.0],
                                 [0.0, self.detector["blue_gain"]]]),
                       (dim_out_2[0], dim_out_2[1]))

        img_raw[0:img.shape[0], 0:img.shape[1]] = (red * img[:, :, 0] +
                                               green * img[:, :, 1] +
                                               blue * img[:, :, 2])

        return img_raw

    def get_demosaiced(self, img_raw, alg="menon"):
        img = DEMOSAIC_ALGOS[alg](img_raw, pattern='RGGB')
        return img

if __name__ == "__main__":
    from matplotlib import pyplot as plt

    d = Detector()

    img = get_snellen_raster(ifov = 900e-6, oversample=2)
    plt.figure()
    plt.imshow(img)

    img_raw = d.get_raw_img(img)
    plt.figure()
    plt.imshow(img_raw)

    img_raw = d.apply_noise_raw(img_raw)
    plt.figure()
    plt.imshow(img_raw)

    img2 = d.get_demosaiced(img_raw)
    d.invert_filter_gain_rgb(img2)
    plt.figure()
    plt.imshow(img2)

    plt.show()
