#!/usr/bin/env python

'''
Make some MUSE stellar & emission line kinematics maps, etc.
for Merry Powell's 2018 paper on the Chandra+MUSE observations of HE0227 and HE0351

G. Tremblay, March 2018
grant.tremblay@cfa.harvard.edu
'''

from astropy.io import fits
import numpy as np

import matplotlib.pyplot as plt

import astropy.constants as const


def main():
	pass

def make_electron_density_map(sii_6716_image, sii_6730_image, hdr):

    ratio = sii_6716_image / sii_6730_image

    # assuming T = 10^4 K, following eq. (3) here: https://arxiv.org/pdf/1311.5041.pdf

    log_ne_per_cm3 = 0.053 * np.tan(-3.0553 * ratio + 2.8506) + \
        6.98 - 10.6905 * ratio + \
        9.9186 * ratio**2 - 3.5442 * ratio**3

    # Mask unrealistic values
    log_ne_mask1 = log_ne_per_cm3 < 0
    log_ne_mask2 = log_ne_per_cm3 > 6

    log_ne_per_cm3[log_ne_mask1] = np.nan
    log_ne_per_cm3[log_ne_mask2] = np.nan

    hdu = fits.PrimaryHDU(log_ne_per_cm3, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('electron_density.fits', overwrite=True,
                    output_verify='silentfix')


def make_balmer_map(ha_image, hb_image, hdr):

    ratio_observed = ha_image / hb_image
    ratio_intrinsic = 2.86
    k_alpha = 2.63
    k_beta = 3.71

    ebv = (2.5 / (k_beta - k_alpha)) * \
        np.log10(ratio_observed / ratio_intrinsic)
    ebv[ebv < 0] = np.nan  # Additional masking

    av = 4.05 * ebv
    nh = 1.8e21 * av  # VERY rough, from Predehl & Schmitt, in atoms / cm2

    hdu = fits.PrimaryHDU(ebv, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('eb-v.fits', overwrite=True, output_verify='silentfix')

    hdu = fits.PrimaryHDU(av, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('a_v.fits', overwrite=True, output_verify='silentfix')

    hdu = fits.PrimaryHDU(nh, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('n_h.fits', overwrite=True, output_verify='silentfix')

    return None


if __name__ == '__main__':
    main()
