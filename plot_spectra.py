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
    
'''
    OneOneFile = '{0}/{0}_NH3_11_{2}_trim.fits'.format(region,file_extension)
    TwoTwoFile = '{0}/{0}_NH3_22_{2}_trim.fits'.format(region,file_extension)
    ThrThrFile = '{0}/{0}_NH3_33_{2}_trim.fits'.format(region,file_extension)



hdu11, header11 = fits.getdata(OneOneFile,header=True)
hdu22, header22 = fits.getdata(TwoTwoFile,header=True)
hdu33, header33 = fits.getdata(ThrThrFile,header=True)
w = WCS(OneOneFile)

# Get spectra and fits
c = SkyCoord('18h30m3.2s','-2d3m6s',frame='icrs')
c11x, c11y = wcsu.skycoord_to_pixel(c,w)
sf = SkyCoord('18h30m12.8s','-2d6m37s',frame='icrs')
sf11x, sf11y = wcsu.skycoord_to_pixel(sf,w)
nc = SkyCoord('18h29m46.7s','-2d6m55s',frame='icrs')
nc11x, nc11y = wcsu.skycoord_to_pixel(nc,w)

cSpec1 = hdu11[:,np.int(c11y),np.int(c11x)]
cFit11 = fit11[:,np.int(c11y),np.int(c11x)]
cSpec2 = hdu22[:,np.int(c11y),np.int(c11x)]
cFit22 = fit22[:,np.int(c11y),np.int(c11x)]
cSpec3 = hdu33[:,np.int(c33y),np.int(c33x)]

sfSpec1 = hdu11[:,np.int(sf11y),np.int(sf11x)]
sfFit11 = fit11[:,np.int(sf11y),np.int(sf11x)]
sfSpec2 = hdu22[:,np.int(sf11y),np.int(sf11x)]
sfFit22 = fit22[:,np.int(sf11y),np.int(sf11x)]
sfSpec3 = hdu33[:,np.int(sf33y),np.int(sf33x)]

ncSpec1 = hdu11[:,np.int(nc11y),np.int(nc11x)]
ncFit11 = fit11[:,np.int(nc11y),np.int(nc11x)]
ncSpec2 = hdu22[:,np.int(nc11y),np.int(nc11x)]
ncFit22 = fit22[:,np.int(nc11y),np.int(nc11x)]
ncSpec3 = hdu33[:,np.int(nc33y),np.int(nc33x)]

# Create velocity axis
# 11 and 22 data:
nchan = header11['NAXIS3']
v_inc = header11['CDELT3']
v_ctr = header11['CRVAL3']
v_pix = header11['CRPIX3']
velo  = (np.arange(nchan)*v_inc + v_inc + (v_ctr - v_inc*v_pix))/1.e3
# 33 data:
nchan = header33['NAXIS3']
v_inc = header33['CDELT3']
v_ctr = header33['CRVAL3']
v_pix = header33['CRPIX3']
velo33 = (np.arange(nchan)*v_inc + v_inc + (v_ctr - v_inc*v_pix))/1.e3

lw = 1.2
fig = plt.figure(figsize=(4,6))
fig.subplots_adjust(hspace=0.1)
ax1 = plt.subplot(3,1,1)
plt.plot(velo,cSpec1,color='black',linewidth=lw,
         label=r'$\mathrm{NH}_3 \ \mathrm{(1,1)}$')
plt.plot(velo,cSpec2-2,color='royalblue',linewidth=lw,
         label=r'$\mathrm{NH}_3 \ \mathrm{(2,2)}$')
plt.plot(velo33,cSpec3*2.-4,color='gray',linewidth=lw,zorder=4,
         label=r'$\mathrm{NH}_3 \ \mathrm{(3,3)} \times \ 2$')
plt.xlim(-18,33)
plt.ylim(-5,10)
plt.yticks(np.arange(-4, 12, 4.0))
ax1.axes.xaxis.set_ticklabels([])
ax1.minorticks_on()
ax1.text(-16,8,'SSC',fontsize=10)
plt.legend(frameon=False,fontsize=8)
ax2 = plt.subplot(3,1,2)
plt.plot(velo,sfSpec1,color='black',linewidth=lw)
plt.plot(velo,sfSpec2-2,color='royalblue',linewidth=lw)
#plt.plot(velo33,sfSpec3,color='gray',linewidth=lw)
plt.xlim(-18,33)
plt.ylim(-3,10)
plt.ylabel(r'$T_\mathrm{MB} \ \mathrm{(K)}$')
plt.yticks(np.arange(0, 12, 4.0))
ax2.axes.xaxis.set_ticklabels([])
ax2.minorticks_on()
ax2.text(-16,8,'S. Filament',fontsize=10)
ax3 = plt.subplot(3,1,3)
plt.plot(velo,ncSpec1,color='black',linewidth=lw)
plt.plot(velo,ncSpec2-1,color='royalblue',linewidth=lw)
#plt.plot(velo33,ncSpec3,color='gray',linewidth=lw)
plt.xlim(-18,33)
plt.ylim(-1.5,3)
plt.yticks(np.arange(-1, 4, 1.0))
plt.xlabel(r'$\mathrm{km \ s}^{-1}$')
ax3.minorticks_on()
ax3.text(-16,2.3,'Narrow-line clump',fontsize=10)
fig.savefig('spectra.pdf',bbox_inches='tight')
plt.close('fig')
'''
