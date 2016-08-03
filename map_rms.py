import matplotlib.pyplot as plt
import astropy.units as u
import warnings
import numpy as np
import os

from astropy.io import fits
import aplpy

from config import plottingDictionary

region_list=['L1688', 'B18', 'NGC1333', 'OrionA']
line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
label_list=['NH$_3$(1,1)','NH$_3$(2,2)','NH$_3$(3,3)','C$_2$S','HC$_5$N',
            'HC$_7$N (21-20)','HC$_7$N (22-21)']
extension='DR1_rebase3'
color_table='Blues'
text_color='black'
beam_color='#d95f02'  # previously used '#E31A1C'
color_hist='#a6bddb'

color_hall= ['#a6cee3', '#fdbf6f', '#33a02c', '#b2df8a', '#fb9a99', '#e31a1c', 
             '#1f78b4']
# color_hall= ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02',
#              '#a6761d']

bin_size=0.01 # K

for region_i in region_list:
    plot_param=plottingDictionary[region_i]#setup_plot_parameters(region=region_i)
    fig_hall=plt.figure(figsize=(6,6))
    ax_hall = fig_hall.add_subplot(111)
    v_hall_max=np.max(plot_param['rms_max'])
    v_hall_min=np.min(plot_param['rms_min'])
    nbin_hall=int(np.ceil( np.abs(v_hall_max-v_hall_min)/ bin_size))
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
            ax.set_xlabel('rms (K)')
            ax.set_ylabel('Fraction of Pixels')
            ax.text(0.9, 0.9, '{0}\n{1}'.format(region_i,label_i),
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes)
            ax.set_xlim(xmin=v_min,xmax=v_max)
            fig.savefig('figures/{0}_{1}_{2}_rms_hist.pdf'.format(region_i,line_i,extension),adjust_bbox=True)
            # now do the histogram of all lines together
            hist, binEdges = np.histogram( myarray, weights=weights, bins=nbin_hall, 
                                       range=(v_hall_min,v_hall_max))
            binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
            ax_hall.plot( binCenters, hist, drawstyle='steps-mid',label=label_i,color=color_hall[i])
        else:
            print('File {0} not found'.format(file_rms))
    # after for loop is finished then close things
    ax_hall.set_xlabel('rms (K)')
    ax_hall.set_ylabel('Fraction of Pixels')
    ax_hall.set_xlim(xmin=v_hall_min,xmax=v_hall_max)
    leg=ax_hall.legend(frameon=False, loc=1, title=region_i)
    leg_text=leg.get_texts()
    for index in range(len(leg_text)):
        plt.setp( leg_text[index], color=color_hall[index] )
    fig_hall.savefig('figures/{0}_all_{1}_rms_hist.pdf'.format(region_i,extension),adjust_bbox=True)
    plt.close('all')
