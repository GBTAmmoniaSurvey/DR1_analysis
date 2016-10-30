import numpy as np
import astropy.io.fits as fits
from spectral_cube import SpectralCube
import astropy.units as u

def deblend_cube(region='OrionA',vmin=0.,vmax=20.):
  tex = fits.getdata('{0}/parameterMaps/{0}_Tex_DR1_rebase3_flag.fits'.format(region))
  mom0 = fits.getdata('{0}/{0}_NH3_11_DR1_rebase3_mom0_QA_trim.fits'.format(region))
  vlsr = fits.getdata('{0}/parameterMaps/{0}_Vlsr_DR1_rebase3_flag.fits'.format(region))
  sigv = fits.getdata('{0}/parameterMaps/{0}_Sigma_DR1_rebase3_flag.fits'.format(region))
  nnh3 = fits.getdata('{0}/parameterMaps/{0}_N_NH3_DR1_rebase3_flag.fits'.format(region))
  cube = SpectralCube.read('{0}/{0}_NH3_11_DR1_rebase3_trim.fits'.format(region))
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
  slab = newcube.spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s)
  slab.write('{0}/{0}_singlecomp.fits'.format(region),overwrite=True)  
      
region_list = ['OrionA','NGC1333','L1688','B18']
vmin_list   = [0.,0.,-4.,0.]
vmax_list   = [20.,12.,10.,12.]

for region_i in range(len(region_list)):
    region = region_list[region_i]
    vmin = vmin_list[region_i]
    vmax = vmax_list[region_i]
    deblend_cube(region,vmin,vmax)




