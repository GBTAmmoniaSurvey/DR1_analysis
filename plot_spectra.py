import astropy.io.fits as fits
from astropy.io.fits import update
from astropy.wcs import WCS
import astropy.wcs.utils as wcsu
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
from astropy.coordinates import SkyCoord
from spectral_cube import SpectralCube

from config_spectra import plottingDictionary
'''
Script to plot some spectra for each region for DR1 paper
For each line show summed or averaged spectrum over entire region
Pick out some interesting locations as well?
'''

def sum_spectra(region,line,file_extension,restFreq,spec_param):
    '''
    Sum spectra over entire mapped region
    Cubes are missing BUNIT header parameter. Fix. 
    '''
    filein = '{0}/{0}_{1}_{2}_trim.fits'.format(region,line,file_extension)
    #add_fits_units(filein,'K')
    cube = SpectralCube.read(filein)
    slice_unmasked = cube.unmasked_data[:,:,:]
    if line == 'NH3_33':
        slice_unmasked[spec_param['mask33_chans'][0]:spec_param['mask33_chans'][1],:,:]=0.
    summed_spectrum = np.nanmean(slice_unmasked,axis=(1,2))
    cube2 = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio',
                                    rest_value=restFreq*u.GHz)
    return summed_spectrum, cube2.spectral_axis

def add_fits_units(filein,bunit):
    hdu = fits.open(filein)
    header = hdu[0].header
    header.set('BUNIT', 'K')
    hdu.writeto(filein,clobber=True)
    hdu.close()

def mask_spikes(cube,tmax):
    mask = cube > tmax * u.K
    cube2 = cube.with_mask(mask)
    return cube2

region_list=['L1688', 'B18', 'NGC1333', 'OrionA']
extension='DR1_rebase3'
# HC7N 22-21 is very noisy
#line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
#label_list=['NH$_3$(1,1)','NH$_3$(2,2)','NH$_3$(3,3)','C$_2$S','HC$_5$N',
#            'HC$_7$N (21-20)','HC$_7$N (22-21)']
line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20']
label_list=['NH$_3$(1,1)','NH$_3$(2,2)','NH$_3$(3,3)','C$_2$S','HC$_5$N',
            'HC$_7$N']
restFreq_list = [23.6944955,23.7226336,23.8701296,22.34403,23.9639010,23.6878974]
xtext_list = [22,22,22,28,28,28]
colour_list = ['black','darkslateblue','blue','royalblue','darkcyan','forestgreen']

for region in region_list:
    spec_param = plottingDictionary[region]
    fig = plt.figure(figsize=(4,7))
    ax = plt.gca()
    plt.xlim(spec_param['vmin'],spec_param['vmax'])
    plt.ylim(spec_param['min_amp'],spec_param['max_amp'])
    ax.set_xlabel('km s$^{-1}$')
    ax.set_ylabel(r'$T_\mathrm{MB}$ (K)')
    y_offset_array = spec_param['y_offsets']
    label_offset_array = spec_param['label_offsets']
    for line_i in range(len(line_list)):
        line = line_list[line_i]
        label = label_list[line_i]
        y_offset = y_offset_array[line_i]
        label_offset = label_offset_array[line_i]
        restFreq = restFreq_list[line_i]
        summed_spectrum,velocity_axis = sum_spectra(region,line,extension,restFreq,spec_param)
        plt.plot(velocity_axis,summed_spectrum.value+y_offset,label=label,color=colour_list[line_i])
        plt.plot([-50,50],[y_offset,y_offset],color='gray',linestyle='dotted',
                 zorder=1,linewidth=0.5)
        
        ax.text(xtext_list[line_i],y_offset+label_offset,label)
    #plt.legend(frameon=False)
    ax.text(0.06,0.95,region,size=16,transform=ax.transAxes)
    ax.minorticks_on()
    plt.tight_layout()
    fig.savefig('figures/{0}_mean_spectra.pdf'.format(region))
    plt.close()

