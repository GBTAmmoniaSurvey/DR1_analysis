import matplotlib.pyplot as plt
import astropy.units as u
import warnings
import numpy as np
import os

from astropy.io import fits
import aplpy

from config import plottingDictionary

"""
Make property map plots for DR1 regions 
Addresses issue #3 in DR1_analysis
"""

#region_list=['L1688', 'B18', 'NGC1333', 'OrionA']
region_list = ['OrionA']
extension='DR1_rebase3'
# Property maps to plot
par_list = ['Vlsr','Sigma','Tkin','Tex','N_NH3']
label_list=['$v_{LSR}$ (km s$^{-1}$)','$\sigma_v$ (km s$^{-1}$)','$T_K$ (K)','$T_{ex}$ (K)','log N(para-NH$_3$)']
# Consider putting these in plotting dictionary for each region
pmin_list = [-5,0,5,2.73,13]
pmax_list = [20,2,40,15,15.5]
# Colour tables for plots
ctable_list = ['YlGnBu_r','YlGnBu_r','hot','hot','hot']
# Contour parameters (currently NH3 moment 0)
cont_color='black'
cont_lw   = 0.6
# Other parameters
text_color='black'
beam_face = 'None'
beam_edge = 'black'
nan_color = '0.95'
#beam_color='#d95f02'  # previously used '#E31A1C'

for region_i in region_list:
    plot_param=plottingDictionary[region_i]
    file_w11='{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region_i,extension)
    cont_levs=2**np.arange( 0,10)*plot_param['w11_step']*4.
    for i in range(len(par_list)):
        par_i = par_list[i]
        label_i = label_list[i]
        par_file = '{0}/parameterMaps/{0}_{1}_{2}_flag.fits'.format(region_i,par_i,extension)
        par_hdu = fits.open(par_file)
        v_min = np.max([pmin_list[i],np.nanmin(par_hdu[0].data)])
        v_max = np.min([pmax_list[i],np.nanmax(par_hdu[0].data)])
        if os.path.isfile(par_file):
            fig=aplpy.FITSFigure(par_file,figsize=(plot_param['size_x'], plot_param['size_y']))
            fig.show_colorscale(cmap=ctable_list[i],vmin=v_min,vmax=v_max)
            fig.set_nan_color(nan_color)
            #
            fig.show_contour(file_w11,colors=cont_color,levels=cont_levs,
                             linewidths=cont_lw)
            #
            # Ticks
            fig.ticks.set_color(text_color)
            fig.tick_labels.set_xformat('hh:mm:ss')
            fig.tick_labels.set_style('colons')
            fig.tick_labels.set_yformat('dd:mm')
            # Add beam
            fig.add_beam()
            fig.beam.set_edgecolor(beam_edge)
            fig.beam.set_facecolor(beam_face)
            fig.beam.set_corner(plot_param['beam_pos'])
            # Scale bar
            # magic line of code to obtain scale in arcsec obtained from 
            # http://www.astropy.org/astropy-tutorials/Quantities.html
            ang_sep = (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
            fig.add_scalebar(ang_sep.to(u.degree))
            fig.scalebar.set_corner(plot_param['scalebar_pos'])
            fig.scalebar.set(color=text_color)
            fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
            # Labels
            fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                          '{0}\n{1}'.format(region_i,label_i), 
                          relative=True, color=text_color, 
                          horizontalalignment=plot_param['label_align'])
            # add colorbar. Colorbar labels not always ideal. 
            fig.add_colorbar()
            fig.colorbar.set_width(0.15)
            fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, location='top')
            # fig.set_system_latex(True)
            fig.save( 'figures/{0}_{1}_{2}_map.pdf'.format(region_i,extension,par_i),adjust_bbox=True)#, bbox_inches='tight')
            fig.close()
        else:
            print('File {0} not found'.format(par_file))

        
        
