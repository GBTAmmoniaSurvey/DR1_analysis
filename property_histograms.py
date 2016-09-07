from astropy.io import fits
import aplpy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import astropy.units as u
import astropy.constants as c
import warnings
import numpy as np
from astropy.visualization import hist
from config import plottingDictionary
"""
Make histogram plots of NH3-derived properties for DR1 regions
"""
def mask_hist(par_data,epar_data,epar_lim,par_max,par_min=0):
    mask1 = np.isfinite(epar_data)
    mask2 = epar_data < epar_lim
    mask3 = epar_data > 0
    mask4 = par_data < par_max
    mask5 = par_data > par_min
    return par_data * mask1 * mask2 * mask3 * mask4 * mask5  


#region_list = ['L1688','NGC1333','B18','OrionA']
region_list = ['B18','NGC1333','L1688','OrionA']
par_list = ['Vlsr','Sigma','Tkin','Tex','N_NH3']
epar_ext = [10,9,6,7,8]
epar_limits = [0.05,0.1,1,2,0.5]
extension = 'DR1_rebase3'
label_list=['$v_{LSR}$ (km s$^{-1}$)','$\sigma_v$ (km s$^{-1}$)','$T_K$ (K)','$T_{ex}$ (K)','log N(para-NH$_3$)']
plot_colours = ['black','blue','green','orange']
#plot_colours = ['black','darkblue','blue','cornflowerblue']
#plot_colours = ['#a6cee3', '#fdbf6f', '#33a02c', '#fb9a99']
hist_minx_list = [2,0,5,2.7,13]
hist_maxx_list = [13,1.5,37,12,15.7]
hist_maxy_list = [2.1,8,0.65,1.6,2]
#hist_maxy_list = [600,800,220,1.5e3,500]
ytick_int_maj = [300,200,50,200,100]
ytick_int_min = [50,50,10,50,20]
dataDir = '/media/DATAPART/projects/GAS/data/'

hist_kwds1 = dict(histtype='stepfilled',alpha=0.2,normed=True)
# All together

fig = plt.figure(figsize=(6.5,8))
for par_i in range(len(par_list)):
    par = par_list[par_i]
    label = label_list[par_i]
    #fig = plt.figure()
    #ax = plt.gca()
    ax = plt.subplot(3,2,par_i+1)
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        plot_param=plottingDictionary[region]
        par_file = dataDir + '{0}/parameterMaps/{0}_{1}_{2}_flag.fits'.format(region,par,extension)
        epar_file = dataDir + '{0}/{0}_parameter_maps_{1}_trim.fits'.format(region,extension)
        epar_hdu = fits.open(epar_file)
        epar_data = epar_hdu[0].data[epar_ext[par_i],:,:]
        epar_hdu.close()
        par_hdu = fits.open(par_file)
        par_data = par_hdu[0].data
        par_hdu.close()
        pmin_list = plot_param['pmin_list']
        pmax_list = plot_param['pmax_list']
        pmin = np.max([pmin_list[par_i],np.nanmin(par_data)])
        pmax = np.min([pmax_list[par_i],np.nanmax(par_data)])
        par_masked = mask_hist(par_data,epar_data,epar_limits[par_i],pmax,par_min=pmin)
        par_masked = par_masked[np.isfinite(par_masked)]
        if par == 'Vlsr':
            bin_width = 0.3
            nbins = np.int((np.max(par_masked) - np.min(par_masked !=0))/bin_width)
            hist(par_masked[par_masked !=0],bins=nbins,ax=ax,histtype='stepfilled',alpha=0.3,
                 normed=True,color=plot_colours[region_i],label=region)
        else:
            hist(par_masked[par_masked !=0],bins='knuth',ax=ax,histtype='stepfilled',alpha=0.3,
                 normed=True,color=plot_colours[region_i],label=region)

    ax.set_xlabel(label)
    ax.set_ylabel('P(t)')
    #ax.set_ylabel('N')
    ax.set_xlim(hist_minx_list[par_i],hist_maxx_list[par_i])
    ax.set_ylim(0,hist_maxy_list[par_i])
    #ax.legend(frameon=False)
    #fig.savefig('figures/{0}_histogram_number.pdf'.format(par))
ax.legend(frameon=False,bbox_to_anchor=(2.0,1.0))
fig.tight_layout()
fig.savefig('figures/all_histograms.pdf')
plt.close('all')

# Separate regions in plots to better show distributions
for par_i in range(len(par_list)):
    fig,axes = plt.subplots(len(region_list),1,figsize=(4.5,5))
    par = par_list[par_i]
    label = label_list[par_i]
    for i, ax in enumerate(fig.axes):
        region_i = i
        region = region_list[region_i]
        plot_param=plottingDictionary[region]
        par_file = dataDir + '{0}/parameterMaps/{0}_{1}_{2}_flag.fits'.format(region,par,extension)
        epar_file = dataDir + '{0}/{0}_parameter_maps_{1}_trim.fits'.format(region,extension)
        epar_hdu = fits.open(epar_file)
        epar_data = epar_hdu[0].data[epar_ext[par_i],:,:]
        epar_hdu.close()
        par_hdu = fits.open(par_file)
        par_data = par_hdu[0].data
        par_hdu.close()
        pmin_list = plot_param['pmin_list']
        pmax_list = plot_param['pmax_list']
        pmin = np.max([pmin_list[par_i],np.nanmin(par_data)])
        pmax = np.min([pmax_list[par_i],np.nanmax(par_data)])
        par_masked = mask_hist(par_data,epar_data,epar_limits[par_i],pmax,par_min=pmin)
        par_masked = par_masked[np.isfinite(par_masked)]
        if par == 'Vlsr':
            bin_width = 0.3
            nbins = np.int((np.max(par_masked) - np.min(par_masked !=0))/bin_width)
            hist(par_masked[par_masked !=0],bins=nbins,ax=ax,histtype='stepfilled',alpha=0.3,
                 color=plot_colours[region_i],label=region)
        else:
            hist(par_masked[par_masked !=0],bins='knuth',ax=ax,histtype='stepfilled',alpha=0.3,
                 color=plot_colours[region_i],label=region)
        if (i+1) != len(region_list):
            ax.set_xticklabels([])
        ax.set_xlim(hist_minx_list[par_i],hist_maxx_list[par_i])
        #ax.set_ylim(0,hist_maxy_list[par_i])
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_int_maj[par_i]))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(ytick_int_min[par_i]))
        ax.annotate('{0}'.format(region),xy=(0.97,0.7),xycoords='axes fraction',horizontalalignment='right')
    #ax.legend(frameon=False)
    ax.set_xlabel(label)
    fig.text(0.01,0.5,'N',va='center',rotation='vertical')
    #fig.tight_layout()
    fig.savefig('figures/{0}_histogram_separated.pdf'.format(par))
    plt.close('all')
