import numpy as np
import astropy.io.fits as fits
from spectral_cube import SpectralCube
import astropy.units as u

region = 'OrionA'

tex = fits.getdata(region+'_Tex_DR1_rebase3_flag.fits')
mom0 = fits.getdata(region+'_NH3_11_DR1_rebase3_mom0_QA_trim.fits')
vlsr = fits.getdata(region+'_Vlsr_DR1_rebase3_flag.fits')
sigv = fits.getdata(region+'_Sigma_DR1_rebase3_flag.fits')
nnh3 = fits.getdata(region+'_N_NH3_DR1_rebase3_flag.fits')
cube = SpectralCube.read(region+'_NH3_11_DR1_rebase3_trim.fits')
cube  = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')

tpeak = mom0/(np.sqrt(2*np.pi)*sigv)

vlsr[vlsr==0]=np.nan
sigv[sigv==0]=np.nan

deblend = np.zeros(cube.shape)
hdr = cube.wcs.to_header()
spaxis = cube.spectral_axis.value

for plane in np.arange(deblend.shape[0]):
    deblend[plane,:,:] = tpeak*np.exp(-(spaxis[plane]-vlsr)**2/(2*sigv**2))

newcube = SpectralCube(deblend,cube.wcs,header=hdr)
newcube.write(region+'_singlecomp.fits',overwrite=True)
