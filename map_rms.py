import matplotlib.pyplot as plt
import astropy.units as u
import warnings
import numpy as np
import os

from astropy.io import fits
import aplpy

region_list=['L1688', 'B18']
line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
label_list=['NH$_3$(1,1)','NH$_3$(2,2)','NH$_3$(3,3)','C$_2$S','HC$_5$N',
            'HC$_7$N (21-20)','HC$_7$N (22-21)']
extension='DR1_rebase3'
color_table='Blues'
text_color='black'
beam_color='#d95f02'  # previously used '#E31A1C'
color_hist='#a6bddb'

bin_size=0.01 # K

def setup_plot_parameters(region='L1688'):
    """
    This function returns a dictionary with the plotting parameters for the 
    requested region.

    There might be a better or more elegant solution, to be discussed

    the rms_max array is the array with the maximum rms to be plotted, where
    the array loops over the lines in hte order of line_list
    """
    if region=='L1688':
        param={'size_x' : 9, 'size_y' : 7, 'scale_bar' : 0.1*u.pc, 
               'distance' : 145.*u.pc, 'beam_pos' : 'bottom left',
               'label_xpos' : 0.025, 'label_ypos' : 0.1,
               'label_align' : 'left',
               'rms_min' : np.array([0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05]),
               'rms_max' : np.array([0.75, 0.75, 0.75, 1.0, 0.75, 0.75, 0.75]) }
    elif region=='B18':
        param={'size_x' : 19, 'size_y' : 7, 'scale_bar' : 0.1*u.pc, 
               'distance' : 145.*u.pc, 'beam_pos' : 'top left',
               'label_xpos' : 0.025, 'label_ypos' : 0.9,
               'label_align' : 'left',
               'rms_min' : np.array([0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05]),
               'rms_max' : np.array([0.75, 0.75, 0.75, 1.0, 0.75, 0.75, 0.75]) }
    else:
        warnings.warn("Region not in DR1 or not defined yet")
    return param


for region_i in region_list:
    plot_param=setup_plot_parameters(region=region_i)
    for i in range(len(line_list)):
        line_i=line_list[i]
        label_i=label_list[i]
        file_rms='{0}/{0}_{1}_{2}_rms.fits'.format(region_i,line_i,extension)
        v_min=plot_param['rms_min'][i]
        v_max=plot_param['rms_max'][i]
        if os.path.isfile(file_rms):
            fig=aplpy.FITSFigure(file_rms, hdu=0, figsize=(plot_param['size_x'], plot_param['size_y']) )
            fig.show_colorscale( cmap=color_table,vmin=v_min, vmax=v_max)
            fig.set_nan_color('0.9')
            #
            # Ticks
            fig.ticks.set_color(text_color)
            fig.tick_labels.set_xformat('hh:mm:ss')
            # Add beam
            fig.add_beam()
            fig.beam.set_color(beam_color)
            fig.beam.set_corner(plot_param['beam_pos'])
            # Scale bar
            # magic line of code to obtain scale in arcsec obtained from 
            # http://www.astropy.org/astropy-tutorials/Quantities.html
            ang_sep =  (plot_param['scale_bar'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
            fig.add_scalebar(ang_sep.to(u.degree))
            fig.scalebar.set(color=text_color)
            fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scale_bar']))
            # Labels
            fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                          '{0}\n{1}'.format(region_i,label_i), 
                          relative=True, color=text_color, 
                          horizontalalignment=plot_param['label_align'])
            # add colorbar
            fig.add_colorbar()
            fig.colorbar.set_width(0.15)
            fig.colorbar.show( box_orientation='horizontal', width=0.1, pad=0.0, 
                                location='top', axis_label_text='(K)')
            # fig.set_system_latex(True)
            fig.save( 'figures/{0}_{1}_{2}_rms_map.pdf'.format(region_i,line_i,extension),adjust_bbox=True)#, bbox_inches='tight')
            fig.close()
            
            data, hd=fits.getdata(file_rms, header=True)
            fig=plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111)
            # the histogram of the data
            nbin=int(np.ceil( np.abs(v_max-v_min)/ bin_size))
            myarray=data[np.isfinite(data)]
            weights = np.ones_like(myarray)/float(len(myarray))
            n, bins, patches = ax.hist( myarray, weights=weights, bins=nbin, 
                                       facecolor=color_hist, range=(v_min,v_max), 
                                       histtype='stepfilled') 
            plt.xlabel('rms (K)')
            plt.ylabel('Fraction of Pixels')
            ax.text(0.9, 0.9, '{0}\n{1}'.format(region_i,label_i),
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes)
            plt.xlim(xmin=v_min,xmax=v_max)
            fig.savefig('figures/{0}_{1}_{2}_rms_hist.pdf'.format(region_i,line_i,extension),adjust_bbox=True)
        else:
            print('File {0} not found'.format(file_rms))

