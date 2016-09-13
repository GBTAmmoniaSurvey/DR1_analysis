import matplotlib.pyplot as plt
import astropy.units as u
import warnings
import numpy as np
import os
from scipy.ndimage import binary_opening

from astropy.io import fits
import aplpy

from config import plottingDictionary

def mask_binary(imageHDU,LowestContour,selem):
    map = imageHDU[0].data
    mask = binary_opening(map > LowestContour, selem)
    MaskedMap = mask*map
    imageHDU[0].data = MaskedMap
    return imageHDU, mask

"""
Make property map plots for DR1 regions 
Addresses issue #3 in DR1_analysis
"""

region_list=['L1688', 'B18', 'NGC1333', 'OrionA']
#region_list = ['NGC1333']
extension='DR1_rebase3'
# Property maps to plot
par_list = ['Vlsr','Sigma','Tkin','Tex','N_NH3']
label_list=['$v_{LSR}$ (km s$^{-1}$)','$\sigma_v$ (km s$^{-1}$)','$T_K$ (K)','$T_{ex}$ (K)','log N(para-NH$_3$)']
# Colour tables for plots
ctable_list = ['RdYlBu_r','YlGnBu_r','hot','hot','hot']
# Contour parameters (currently NH3 moment 0)
cont_color='black'
cont_lw   = 0.6
# Other parameters
text_color='black'
beam_face = 'None'
beam_edge = 'black'
nan_color = '0.95'
#beam_color='#d95f02'  # previously used '#E31A1C'
# Masking of small (noisy) regions
# Property maps have already been flagged
selem = np.array([[0,1,0],[1,1,1],[0,1,0]])
maskLim = [0.25,0.5,0.5,0.5]

for region_i in range(len(region_list)):
    region = region_list[region_i]
    plot_param=plottingDictionary[region]
    # Removing lonely pixels from the NH3 (1,1) contours 
    file_w11='{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,extension)
    cont_levs=2**np.arange( 0,10)*plot_param['w11_step']
    LowestContour = cont_levs[0]*maskLim[region_i]
    w11_hdu = fits.open(file_w11)
    w11_hdu, mask = mask_binary(w11_hdu,LowestContour,selem)
    pmin_list = plot_param['pmin_list']
    pmax_list = plot_param['pmax_list']
    p_step    = plot_param['p_step']
    for par_i in range(len(par_list)):
        par = par_list[par_i]
        label_i = label_list[par_i]
        par_file = '{0}/parameterMaps/{0}_{1}_{2}_flag.fits'.format(region,par,extension)
        if os.path.isfile(par_file):
            par_hdu = fits.open(par_file)
            par_data = par_hdu[0].data
            par_hdu[0].data[par_hdu[0].data == 0] = np.nan
            v_min = np.max([pmin_list[par_i],np.nanmin(par_hdu[0].data)])
            v_max = np.min([pmax_list[par_i],np.nanmax(par_hdu[0].data)])
            fig=aplpy.FITSFigure(par_hdu,figsize=(plot_param['size_x'], plot_param['size_y']))
            fig.show_colorscale(cmap=ctable_list[par_i],vmin=v_min,vmax=v_max)
            fig.set_nan_color(nan_color)
            # Moment contours
            fig.show_contour(w11_hdu,colors=cont_color,levels=cont_levs,
                             linewidths=cont_lw)
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
                          '{0}\n{1}'.format(region,label_i), 
                          relative=True, color=text_color, 
                          horizontalalignment=plot_param['label_align'])
            # Add colorbar with defined tick labels
            colorbar_step = p_step[par_i]
            # make nice tick labels using defined colorbar step
            ticks = np.arange(np.ceil((v_max-v_min)/colorbar_step)+2)*colorbar_step + np.int(v_min)
            fig.add_colorbar()
            fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, location='top',ticks=ticks)
            # fig.set_system_latex(True)
            fig.save( 'figures/{0}_{1}_{2}_map.pdf'.format(region,extension,par),adjust_bbox=True,dpi=100)#, bbox_inches='tight')
            fig.close()
        else:
            print('File {0} not found'.format(par_file))

        
        
