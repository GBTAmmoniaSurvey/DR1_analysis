import matplotlib.pyplot as plt
import astropy.units as u
import warnings
import numpy as np
import os
from scipy.ndimage import binary_opening

from astropy.io import fits
import aplpy

from config import plottingDictionary

region_list=['L1688', 'B18', 'NGC1333', 'OrionA']
#region_list=['OrionA']
line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
label_list=['NH$_3$(1,1)','NH$_3$(2,2)','NH$_3$(3,3)','C$_2$S','HC$_5$N',
            'HC$_7$N (21-20)','HC$_7$N (22-21)']
extension='DR1_rebase3'
color_table='magma'
text_color='black'
beam_color='#d95f02'  # previously used '#E31A1C'
# Masking of small (noisy) regions
selem = np.array([[0,1,0],[1,1,1],[0,1,0]])

for region_i in region_list:
    plot_param=plottingDictionary[region_i]
    file_w11='{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region_i,extension)
    cont_levs=2**np.arange( 0,20)*plot_param['w11_step']
    LowestContour = cont_levs[0]*0.5
    w11_hdu = fits.open(file_w11)
    map = w11_hdu[0].data
    mask = binary_opening(map > LowestContour, selem)
    MaskedMap = mask*map
    w11_hdu[0].data = MaskedMap
    for i in range(len(line_list)):
        line_i=line_list[i]
        label_i=label_list[i]
        file_mom0='{0}/{0}_{1}_{2}_mom0_QA_trim.fits'.format(region_i,line_i,extension)
        v_min=plot_param['mom0_min'][i]
        v_max=plot_param['mom0_max'][i]
        if os.path.isfile(file_mom0):
            fig=aplpy.FITSFigure(file_mom0, hdu=0, figsize=(plot_param['size_x'], plot_param['size_y']) )
            if line_i == 'NH3_11':
                fig.show_colorscale( cmap=color_table,vmin=v_min, vmax=v_max, stretch=plot_param['mom0_stretch'],vmid=v_min-(1.*np.abs(v_min)))
                cbar_ticks = [0,3,6,12,24,48,96]
                # add colorbar
                fig.add_colorbar()
                #fig.colorbar.set_width(0.15)
                fig.colorbar.show( box_orientation='horizontal', width=0.1, pad=0.0, ticks=cbar_ticks,
                                   location='top', axis_label_text='Integrated Intensity (K km s$^{-1}$)')
            else:
                fig.show_colorscale( cmap=color_table,vmin=v_min, vmax=v_max)
                # add colorbar
                fig.add_colorbar()
                #fig.colorbar.set_width(0.15)
                fig.colorbar.show( box_orientation='horizontal', width=0.1, pad=0.0, 
                                   location='top', axis_label_text='Integrated Intensity (K km s$^{-1}$)')
            fig.set_nan_color('0.95')
            #
            fig.show_contour(w11_hdu, colors='gray', levels=cont_levs)
            #
            # Ticks
            fig.ticks.set_color(text_color)
            fig.tick_labels.set_style('colons')
            fig.tick_labels.set_xformat('hh:mm:ss')
            fig.tick_labels.set_yformat('dd:mm')
            # Add beam
            fig.add_beam(major=0.0088441,minor=0.0088441,angle=0)
            fig.beam.set_color(beam_color)
            fig.beam.set_corner(plot_param['beam_pos'])
            # Scale bar
            # magic line of code to obtain scale in arcsec obtained from 
            # http://www.astropy.org/astropy-tutorials/Quantities.html
            ang_sep =  (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
            fig.add_scalebar(ang_sep.to(u.degree))
            fig.scalebar.set_corner(plot_param['scalebar_pos'])
            fig.scalebar.set(color=text_color)
            fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
            # Labels
            fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                          '{0}\n{1}'.format(region_i,label_i), 
                          relative=True, color=text_color, 
                          horizontalalignment=plot_param['label_align'])
            # fig.set_system_latex(True)
            fig.save( 'figures/{0}_{1}_{2}_mom0_map.pdf'.format(region_i,line_i,extension),adjust_bbox=True)
            fig.close()
        else:
            print('File {0} not found'.format(file_mom0))
