import numpy as np
import aplpy
import matplotlib.pyplot as plt
from astropy.io import fits
from config import plottingDictionary

def mask_obs_area(rms_data):
    rms_data[np.isfinite(rms_data)] = 1.
    rms_data[np.isnan(rms_data)] = 0.
    return rms_data

def make_obs_mask(region_list,file_ext):
    for region in region_list:
        rms_hdu = fits.open('{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_ext))
        # Set nan pixels to zero to create mask
        obs_mask = mask_obs_area(rms_hdu[0].data)
        new_hdu = fits.PrimaryHDU(obs_mask,rms_hdu[0].header)
        new_hdu.writeto('{0}/{0}_NH3_11_{1}_obsMask.fits'.format(region,file_ext),clobber=True)

region_list = ['B18','NGC1333','L1688','OrionA']
file_ext    = 'DR1_rebase3'
make_obs_mask(region_list,file_ext)
