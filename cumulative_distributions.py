from astropy.io import fits
import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
import warnings
import numpy as np
import os
import FITS_tools

'''
Scripts to calculate cumulative distributions of H2, fractional detections of NH3 in DR1 regions
To re-run, need to revise location of NH3 moment maps, H2 column density files in main, regrid_h2
'''

def main():
    region_list  = ['L1688','NGC1333','B18','OrionA']
    herFile_list=['OphL1688','perseus','Tau_B18','orionA-N']
    projDir = '/media/DATAPART/projects/GAS/testing/'
    herDir = '/media/DATAPART/projects/GAS/otherData/herschel_ayushi/'
    file_extension='DR1_rebase3'
    regrid_h2(projDir=projDir,region_list=region_list,file_extension=file_extension,
              herDir = herDir, herFile_list=herFile_list)
    calc_cumulative_dist(momlim=0.45,region_list=region_list,file_extension=file_extension)
    cdf_h2(region_list=region_list)
    cumulative_two_panel_sims(momlim=0.45,region_list=region_list)
    cumulative_three_panel_sims(momlim=0.45,region_list=region_list)


def regrid_h2(projDir='/media/DATAPART/projects/GAS/testing/',
              region_list = ['L1688','NGC1333','B18','OrionA'],
              file_extension='DR1_rebase3',
              herDir = '/media/DATAPART/projects/GAS/otherData/herschel_ayushi/',
              herFile_list=['OphL1688','perseus','Tau_B18','orionA-N']):
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        herFilename = herFile_list[region_i]
        herColFits = herDir+'{0}/Colden_temp/{1}_colden_masked.fits'.format(region,herFilename)
        nh3ImFits  = projDir + '{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension)
        h2_hdu  = fits.open(herColFits)
        nh3_hdr = fits.getheader(nh3ImFits)
        new_h2  = FITS_tools.hcongrid.hcongrid(h2_hdu[0].data,h2_hdu[0].header,nh3_hdr)
        new_h2_hdu = fits.PrimaryHDU(new_h2,nh3_hdr)
        new_h2_hduList = fits.HDUList([new_h2_hdu])
        new_h2_hduList.writeto('nh2_regridded/{0}_NH2_regrid.fits',clobber=True)

def calc_cumulative_dist(momlim=0.3,region_list=['L1688','NGC1333','B18','OrionA'],
                         file_extension='DR1_rebase3'):
    projDir = '/media/DATAPART/projects/GAS/testing/'
    min_bin = 20.5 # log 10 N(H2) - set by B18
    max_bin = 26. # log 10 N(H2) - set by Orion A
    bin_size = 0.05
    nbins = np.int((max_bin - min_bin)/bin_size)
    data_bin_vals = [min_bin + x * bin_size for x in range(nbins+1)]
    # Loop over regions
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        nh3ImFits  = projDir + '{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension)
        herColFits = 'nh2_regridded/{0}_NH2_regrid.fits'.format(region)
        nh3_fits = fits.open(nh3ImFits)
        nh2_regrid_fits = fits.open(herColFits)
        h2_data = np.log10(nh2_regrid_fits[0].data)
        nh3_data = nh3_fits[0].data
        h2_mean_array = np.zeros(nbins)
        h2_std_array = np.zeros(nbins)
        nh3_frac_array = np.zeros(nbins)
        for bin_i in range(nbins-1):
            bin_h2_indices = np.where(np.logical_and(h2_data >= data_bin_vals[bin_i],
                                                     h2_data < data_bin_vals[bin_i+1]))
            bin_h2_data = h2_data[bin_h2_indices]
            bin_nh3_data = nh3_data[bin_h2_indices]
            if np.count_nonzero(bin_nh3_data) != 0:
                frac_above_mom = np.count_nonzero(bin_nh3_data > momlim)/(1.0*np.count_nonzero(bin_nh3_data))
            else:
                frac_above_mom = 0
            bin_mean_h2 = np.nanmean(bin_h2_data)
            bin_std_h2 = np.nanstd(bin_h2_data)
            h2_mean_array[bin_i] = bin_mean_h2
            nh3_frac_array[bin_i] = frac_above_mom
            h2_std_array[bin_i] = bin_std_h2
            # Write out bins
            np.savetxt('cumulative/{0}_cumulative.txt'.format(region),
                       np.transpose([h2_mean_array,nh3_frac_array,h2_std_array]))


def cumulative_two_panel_sims(momlim=0.3,region_list=['L1688','NGC1333','B18','OrionA']):
    # Plot Stella's simulations, include regions as background curves
    # Revise to do double plot here and single plot showing full extent of N(H2) above
    sim_file_total = 'sims/Column_total_45.txt'
    sim_file_hi    = 'sims/Column_hi_45.txt'
    sim_file_lo    = 'sims/Column_lo4_45.txt'
    nh2_tot, fnh3_tot = np.loadtxt(sim_file_total,usecols=(0,1),skiprows=1,unpack=True)
    nh2_hi, fnh3_hi = np.loadtxt(sim_file_hi,usecols=(0,1),skiprows=1,unpack=True)
    nh2_lo, fnh3_lo = np.loadtxt(sim_file_lo,usecols=(0,1),skiprows=1,unpack=True)
# Read in bin data for regions
    b18_file      = 'cumulative/B18_cumulative.txt'
    ngc1333_file  = 'cumulative/NGC1333_cumulative.txt'
    l1688_file    = 'cumulative/L1688_cumulative.txt'
    oriona_file   = 'cumulative/OrionA_cumulative.txt'
    file_list = [l1688_file,ngc1333_file,b18_file,oriona_file]
# Make zoomed in plot only
    min_bin = 20.5 # log 10 N(H2) - set by B18
    max_bin = 23. # log 10 N(H2) - set by Orion A
# Set up plot
    fig = plt.figure()
    ax = plt.subplot(2,1,1)
    ax2 = plt.subplot(2,1,2)
    ax.set_xlim(min_bin, max_bin)
    ax2.set_xlim(min_bin,max_bin)
    ax.set_ylim(-0.05,1.1)
    ax2.set_ylim(-0.05,1.1)
    ax2.set_xlabel('log N(H$_2$) (cm$^{-2}$)')
    plot_colours = ['black','darkblue','blue','cornflowerblue']
    plot_symbols = ['o','v','^','s']
# Loop over regions
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        nh3_filein = file_list[region_i]
        nh2_bin, fnh3_bin = np.loadtxt(nh3_filein,usecols=(0,1),unpack=True)
    # Plot
        ax.scatter(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],
                   marker=plot_symbols[region_i],c=plot_colours[region_i],label=region,s=30)
        ax.plot(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],zorder=1,
                linewidth=2,alpha=0.3,color=plot_colours[region_i],label='__nolabel__')
        ax2.plot(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],
                 linewidth=2,alpha=0.3,color=plot_colours[region_i],label='__nolabel__')
# Add simulation data
    ax2.scatter(nh2_tot,fnh3_tot,marker=plot_symbols[0],c=plot_colours[0],s=30,label=r'$900\ \mathrm{cm}^{-3}$')
    ax2.scatter(nh2_hi,fnh3_hi,marker=plot_symbols[1],c=plot_colours[1],s=30,label=r'$1700\ \mathrm{cm}^{-3}$')
    ax2.scatter(nh2_lo,fnh3_lo,marker=plot_symbols[2],c=plot_colours[2],s=30,label=r'$830\ \mathrm{cm}^{-3}$')
    ax.plot([21.82,21.82],[-1,2],linestyle='--',linewidth=4,color='gray',alpha=0.5,zorder=1)
    ax2.plot([21.82,21.82],[-1,2],linestyle='--',linewidth=4,color='gray',alpha=0.5,zorder=1)
    fig.text(0.015,0.5,'Fraction of pixels above {0} K km/s'.format(momlim),
             va='center',rotation='vertical',fontsize=14)
    fig.text(0.05,0.3,'Simulations',va='center',rotation='vertical')
    fig.text(0.05,0.75,'Observations',va='center',rotation='vertical')
    ax.legend(frameon=False,loc=2,fontsize=12,scatterpoints=1)
    ax2.legend(frameon=False,loc=2,fontsize=14,scatterpoints=1)
    plt.savefig('figures/allDR1_nh3_nh2_cumulative_withSims.pdf')
    plt.close()

def cdf_h2(region_list=['L1688','NGC1333','B18','OrionA']):
    # Calculate the CDF of the N(H2) values for all the DR1 regions
    projDir = '/media/DATAPART/projects/GAS/'
    region_list  = ['L1688','NGC1333','B18','OrionA']
    file_extension = 'DR1_rebase3'
    min_bin = 20.5 # log 10 N(H2) - set by B18
    max_bin = 26. # log 10 N(H2) - set by Orion A
    bin_size = 0.05
    nbins = np.int((max_bin - min_bin)/bin_size)
    data_bin_vals = np.array([min_bin + x * bin_size for x in range(nbins)])
    # Loop over regions
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        # Use regridded files to match NH3 image extent
        herColFits = 'nh2_regridded/{0}_NH2_regrid.fits'.format(region)
        herCol = fits.open(herColFits)
        # Need to mask where no NH3 observations
        nh3ImFits  = projDir + 'data/{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension)
        nh3Im = fits.open(nh3ImFits)
        h2_data = np.log10(herCol[0].data)
        h2_data[np.isnan(nh3Im[0].data)] = 0.
        h2_cumulative_array = np.zeros(nbins)
        for bin_i in range(nbins-1):
            bin_h2_indices = np.where(h2_data < data_bin_vals[bin_i])
            bin_h2_data = h2_data[bin_h2_indices]
            frac_below_bin = np.count_nonzero(bin_h2_data)/(1.0*np.count_nonzero(h2_data))
            h2_cumulative_array[bin_i] = frac_below_bin
        # Write out bins
        np.savetxt('cumulative/{0}_cumulative_nh2.txt'.format(region),
                   np.transpose([data_bin_vals,h2_cumulative_array]))

def cumulative_three_panel_sims(momlim=0.3,region_list=['L1688','NGC1333','B18','OrionA']):
    # Make cumulative distribution function
    cdf_h2()
    # Plot Stella's simulations, include regions as background curves
    # Revise to do triple plot here 
    # Add in CDF of column density PDF
    sim_file_total = 'sims/Column_total_45.txt'
    sim_file_hi    = 'sims/Column_hi_45.txt'
    sim_file_lo    = 'sims/Column_lo4_45.txt'
    nh2_tot, fnh3_tot = np.loadtxt(sim_file_total,usecols=(0,1),skiprows=1,unpack=True)
    nh2_hi, fnh3_hi = np.loadtxt(sim_file_hi,usecols=(0,1),skiprows=1,unpack=True)
    nh2_lo, fnh3_lo = np.loadtxt(sim_file_lo,usecols=(0,1),skiprows=1,unpack=True)
# Make zoomed in plot only
    min_bin = 20.5 # log 10 N(H2) - set by B18
    max_bin = 23. # log 10 N(H2) - set by Orion A
# Set up plot
    fig = plt.figure(figsize=(5,7))
    ax = plt.subplot(3,1,1)
    ax2 = plt.subplot(3,1,2)
    ax3 = plt.subplot(3,1,3)
    ax.set_xlim(min_bin, max_bin)
    ax2.set_xlim(min_bin,max_bin)
    ax3.set_xlim(min_bin,max_bin)
    ax.set_ylim(-0.05,1.1)
    ax2.set_ylim(-0.05,1.1)
    ax3.set_ylim(-0.05,1.1)
    ax3.set_xlabel('log N(H$_2$) (cm$^{-2}$)')
    plot_colours = ['black','darkblue','blue','cornflowerblue']
    plot_symbols = ['o','v','^','s']
# Loop over regions
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        nh3_filein = 'cumulative/{0}_cumulative.txt'.format(region)
        nh2_bin, fnh3_bin = np.loadtxt(nh3_filein,usecols=(0,1),unpack=True)
        # Plot
        ax.scatter(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],
                   marker=plot_symbols[region_i],c=plot_colours[region_i],label=region,s=20)
        ax.plot(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],zorder=1,
                linewidth=2,alpha=0.3,color=plot_colours[region_i],label='__nolabel__')
        ax2.plot(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],
                 linewidth=2,alpha=0.3,color=plot_colours[region_i],label='__nolabel__')
        nh2_filein = 'cumulative/{0}_cumulative_nh2.txt'.format(region)
        nh2_bin,fnh2 = np.loadtxt(nh2_filein,usecols=(0,1),unpack=True)
        # Plot
        ax3.scatter(nh2_bin,fnh2,marker=plot_symbols[region_i],c=plot_colours[region_i],label=region,s=20)
        ax3.plot(nh2_bin,fnh2,zorder=1,linewidth=2,alpha=0.3,color=plot_colours[region_i],label='__nolabel__')
# Add simulation data
    ax2.scatter(nh2_tot,fnh3_tot,marker=plot_symbols[0],c=plot_colours[0],s=30,label=r'$900\ \mathrm{cm}^{-3}$')
    ax2.scatter(nh2_hi,fnh3_hi,marker=plot_symbols[1],c=plot_colours[1],s=30,label=r'$1700\ \mathrm{cm}^{-3}$')
    ax2.scatter(nh2_lo,fnh3_lo,marker=plot_symbols[2],c=plot_colours[2],s=30,label=r'$830\ \mathrm{cm}^{-3}$')
    ax.plot([21.82,21.82],[-1,2],linestyle='--',linewidth=4,color='gray',alpha=0.5,zorder=1)
    ax2.plot([21.82,21.82],[-1,2],linestyle='--',linewidth=4,color='gray',alpha=0.5,zorder=1)
    ax3.plot([21.82,21.82],[-1,2],linestyle='--',linewidth=4,color='gray',alpha=0.5,zorder=1)
    fig.text(0.04,0.7,'Fraction of pixels > {0} K km/s'.format(momlim),
             va='center',rotation='vertical',fontsize=13)
    fig.text(0.65,0.47,'Simulations',va='center',fontsize=11)
    fig.text(0.65,0.8,'Observations',va='center',fontsize=11)
    ax.legend(frameon=False,loc=2,fontsize=11,scatterpoints=1)
    ax2.legend(frameon=False,loc=2,fontsize=13,scatterpoints=1)
    ax3.legend(frameon=False,loc=4,fontsize=11,scatterpoints=1)
    ax3.set_ylabel('Fraction of pixels < N(H$_2$)',fontsize=13)
    plt.tight_layout()
    plt.savefig('figures/allDR1_nh3_nh2_cumulative_withSims_3panel.pdf')
    plt.close()

def cumulative_one_panel(momlim=0.3,region_list=['L1688','NGC1333','B18','OrionA']):
    # New single plot for full N(H2) extent
    min_bin = 20.5 # log 10 N(H2) - set by B18
    max_bin = 26. # log 10 N(H2) - set by Orion A
# Set up plot
    fig = plt.figure()
    ax = plt.subplot(2,1,1)
# Loop over regions
    for region_i in range(len(region_list)):
        region = region_list[region_i]
        nh3_filein = file_list[region_i]
        nh2_bin, fnh3_bin = np.loadtxt(nh3_filein,usecols=(0,1),unpack=True)
    # Plot
        ax.scatter(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],
                   marker=plot_symbols[region_i],c=plot_colours[region_i],label=region,s=30)
        ax.plot(nh2_bin[nh2_bin != 0],fnh3_bin[nh2_bin != 0],zorder=1,
                linewidth=2,alpha=0.3,color=plot_colours[region_i],label='__nolabel__')
    ax.add_patch(patches.Rectangle((np.log10(6.e21),-0.05),0.07,1.2,alpha=0.2,zorder=1))
    ax.set_ylim(-0.05,1.1)
    ax.set_xlim(min_bin,max_bin)
    ax.set_xlabel('log N(H$_2$) (cm$^{-2}$)')
    fig.text(0.015,0.71,'Fraction of pixels above {0} K km/s in NH$_3$'.format(momlim),
             va='center',rotation='vertical',fontsize=14)
    fig.text(0.05,0.71,'SNR=3 in NH$_3$',va='center',rotation='vertical',fontsize=14)
    ax.legend(frameon=False,loc=4,fontsize=12,scatterpoints=1)
    plt.savefig('figures/allDR1_nh3_nh2_cumulative_single.pdf')
    plt.close()

main()
