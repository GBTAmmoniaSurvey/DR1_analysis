from astropy.io import fits
import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
import warnings
import numpy as np
import os
import FITS_tools
import matplotlib.patches as patches
from scipy.ndimage import binary_opening
from GAS_hers_analysis_config import plottingDictionary

def get_prot_loc(region_name):
    # C2D doesn't include Orion, Taurus
    c2d_data_file = 'protostars/c2d_dunham_table2_edit.txt'
    orion_file    = 'protostars/megeath_orion_table4_edit.txt'
    taurus1_file  = 'protostars/rebull_taurus_table6_edit.txt'
    taurus2_file  = 'protostars/rebull_taurus_table7_edit.txt'
    if region_name == 'L1688' or region_name == 'NGC1333':
        protData = np.genfromtxt(c2d_data_file,usecols=(1,2,3,4,5,6,7,12,15),
                                 dtype={'names': ('Region',
                                                  'ra1','ra2','ra3','de1','de2','de3',
                                                  'alpha','AGB'),
                                        'formats': ('S16','i','i','d','i','i','d','d','S1')},
                                 unpack=True,skip_header=24)
        ra_deg = 15. * (protData['ra1']+(protData['ra2'] + protData['ra3']/60.)/60.)
        de_deg = np.abs(protData['de1'])+(protData['de2'] + protData['de3']/60.)/60.
        de_deg[protData['de1'] < 0] = de_deg[protData['de1']<0] * (-1.)
        # For Class 0+I
        #ra_deg_class0I = ra_deg[np.logical_and(protData['AGB'] == 'N',protData['alpha'] >= 0.3)]
        #de_deg_class0I = de_deg[np.logical_and(protData['AGB'] == 'N',protData['alpha'] >= 0.3)]
        # For Class 0+I and flat spectrum
        ra_deg_class0I = ra_deg[np.logical_and(protData['AGB'] == 'N',protData['alpha'] >= -0.3)]
        de_deg_class0I = de_deg[np.logical_and(protData['AGB'] == 'N',protData['alpha'] >= -0.3)]
    if region_name == 'OrionA':
        protData = np.genfromtxt(orion_file,usecols=(1,2,3,4,5,6,7),
                                 dtype={'names':('ra1','ra2','ra3','de1','de2','de3','class'),
                                        'formats':('i','i','i','i','i','i','S5')},
                                 unpack=True,skip_header=51)
        ra_deg = 15. * (protData['ra1']+(protData['ra2'] + protData['ra3']/60.)/60.)
        de_deg = np.abs(protData['de1'])+(protData['de2'] + protData['de3']/60.)/60.
        de_deg[protData['de1'] < 0] = de_deg[protData['de1']<0] * (-1.)
        ra_deg_class0I = ra_deg[protData['class'] == 'P']
        de_deg_class0I = de_deg[protData['class'] == 'P']
    if region_name == 'B18':
        protData = np.genfromtxt(taurus2_file,usecols=(0,1,2,3,4,5,6,7),
                                 dtype={'names':('ra1','ra2','ra3','de1','de2','de3','grade','class'),
                                        'formats':('i','i','i','i','i','i','S2','S5')},
                                 unpack=True,skip_header=71)
        ra_deg = 15. * (protData['ra1']+(protData['ra2'] + protData['ra3']/60.)/60.)
        de_deg = np.abs(protData['de1'])+(protData['de2'] + protData['de3']/60.)/60.
        de_deg[protData['de1'] < 0] = de_deg[protData['de1']<0] * (-1.)
        ra_deg_class0I_1 = ra_deg[np.logical_and(protData['grade'] == 'A',
                                               np.logical_or(protData['class'] == 'I',
                                                             protData['class'] == 'flat'))]
        de_deg_class0I_1 = de_deg[np.logical_and(protData['grade'] == 'A',
                                               np.logical_or(protData['class'] == 'I',
                                                             protData['class'] == 'flat'))]
        protData2 = np.genfromtxt(taurus1_file,usecols=(0,1,2,3,4,5,6),
                                  dtype={'names':('ra1','ra2','ra3','de1','de2','de3','class'),
                                         'formats':('i','i','i','i','i','i','S5')},
                                  unpack=True,skip_header=37)
        ra_deg = 15. * (protData2['ra1']+(protData2['ra2'] + protData2['ra3']/60.)/60.)
        de_deg = np.abs(protData2['de1'])+(protData2['de2'] + protData2['de3']/60.)/60.
        de_deg[protData2['de1'] < 0] = de_deg[protData2['de1']<0] * (-1.)
        ra_deg_class0I_2 = ra_deg[np.logical_or(protData2['class'] == 'I',
                                                protData2['class'] == 'flat')]
        de_deg_class0I_2 = de_deg[np.logical_or(protData2['class'] == 'I',
                                                protData2['class'] == 'flat')]
        ra_deg_class0I = np.append(ra_deg_class0I_1,ra_deg_class0I_2)
        de_deg_class0I = np.append(de_deg_class0I_1,de_deg_class0I_2)
    
    if not os.path.isfile('protostars/{0}_protostar_list.txt'.format(region_name)):
        np.savetxt('protostars/{0}_protostar_list.txt'.format(region_name),
                   np.transpose([ra_deg_class0I,de_deg_class0I]))
    return ra_deg_class0I, de_deg_class0I

def regrid_h2(nh3_image,h2_image):
    # Edit to write out regridded image - glue won't work if files not on same grid
    h2fits  = fits.open(h2_image)
    nh3_hdr = fits.getheader(nh3_image)
    new_h2  = FITS_tools.hcongrid.hcongrid(h2fits[0].data,h2fits[0].header,nh3_hdr)
    new_h2_hdu = fits.PrimaryHDU(new_h2,nh3_hdr)
    return new_h2_hdu

def log10_h2(h2_image):
    log_h2_data = np.log10(h2_image.data)
    log_h2_hdu = fits.PrimaryHDU(log_h2_data,h2_image.header)
    return log_h2_hdu    

def plot_dust_overlay(nh3_image_fits,h2_image_fits,region,plot_param,v_min,v_max,maskLim,obsMaskFits):
    text_size = 14
    b18_text_size = 20   
    # Contour parameters (currently NH3 moment 0)
    cont_color='black'
    cont_lw   = 0.6
    cont_levs=2**np.arange( 0,20)*plot_param['w11_step']
    # Masking of small (noisy) regions
    selem = np.array([[0,1,0],[1,1,1],[0,1,0]])
    LowestContour = cont_levs[0]*0.5
    w11_hdu = fits.open(nh3_image_fits)
    map = w11_hdu[0].data
    mask = binary_opening(map > LowestContour, selem)
    MaskedMap = mask*map
    w11_hdu[0].data = MaskedMap
    # Labels
    if region == 'B18':
        label_colour = 'white'
        text_size = b18_text_size
    else:
        label_colour = 'white'
    fig=aplpy.FITSFigure(h2_image_fits,figsize=(plot_param['size_x'], plot_param['size_y']))
    fig.show_colorscale(cmap='hot',vmin=v_min,vmax=v_max,stretch='log',vmid=v_min-v_min*plot_param['vmid_scale'])
    fig.set_nan_color('0.95')
    # Observations mask contour
    fig.show_contour(obsMaskFits,colors='white',levels=1,linewidths=1.5)
    # NH3 moment contours
    fig.show_contour(w11_hdu,colors=cont_color,levels=cont_levs,linewidths=cont_lw)
    fig.axis_labels.set_font(family='sans_serif',size=text_size)
    # Ticks
    fig.tick_labels.set_font(family='sans_serif',size=text_size)
    fig.ticks.set_color('white')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_style('colons')
    fig.tick_labels.set_yformat('dd:mm')
    # Scale bar
    # magic line of code to obtain scale in arcsec obtained from 
    # http://www.astropy.org/astropy-tutorials/Quantities.html
    ang_sep = (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
    fig.add_scalebar(ang_sep.to(u.degree))
    fig.scalebar.set_font(family='sans_serif',size=text_size)
    fig.scalebar.set_corner(plot_param['scalebar_pos'])
    fig.scalebar.set(color='white')
    fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
    fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                  '{0}\n{1}'.format(region,'500 $\mu$m'), 
                  relative=True, color=label_colour, 
                  horizontalalignment=plot_param['label_align'],
                  family='sans_serif',size=text_size)

    fig.save( 'figures/{0}_continuum_image.pdf'.format(region),adjust_bbox=True,dpi=200)#, bbox_inches='tight')
    fig.close()

def plot_column_overlay(nh3_cont_fits,h2_col_fits,region,plot_pars,maskLim,obsMaskFits,include_prot=True):
    fix_text_size = 14
    b18_text_size = 20    
    # Color scale parameters (for log(data))
    vmin = 20.5
    vmax = 23.
    # Contour parameters (currently NH3 moment 0)
    cont_color='indigo'
    cont_lw   = 0.6
    cont_levs=2**np.arange( 0,20)*plot_param['w11_step']
    fig=aplpy.FITSFigure(h2_col_fits,figsize=(plot_param['size_x'], plot_param['size_y']))
    fig.show_colorscale(cmap='afmhot',vmin=vmin,vmax=vmax)
    fig.set_nan_color('0.95')
    # Add outline of observations 
    fig.show_contour(nh3_cont_fits,colors='white',levels=[np.nan],linewidths=cont_lw)
    # Observations mask contour
    fig.show_contour(obsMaskFits,colors='white',levels=1,linewidths=1.5)
    # NH3 moment contours
    # Masking of small (noisy) regions
    selem = np.array([[0,1,0],[1,1,1],[0,1,0]])
    LowestContour = cont_levs[0]*0.5
    w11_hdu = fits.open(nh3_cont_fits)
    map = w11_hdu[0].data
    mask = binary_opening(map > LowestContour, selem)
    MaskedMap = mask*map
    w11_hdu[0].data = MaskedMap
    # Labels & colours
    if region == 'B18':
        label_colour = 'black'
        text_size = b18_text_size
        tick_colour = 'black'
    elif region == 'NGC1333':
        label_colour = 'white'
        tick_colour = 'white'
        text_size = b18_text_size
    else:
        label_colour = 'white'
        tick_colour = 'white'
        text_size = fix_text_size
    fig.show_contour(w11_hdu,colors=cont_color,levels=cont_levs,linewidths=cont_lw)
    # Ticks
    fig.ticks.set_color(tick_colour)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_style('colons')
    fig.tick_labels.set_yformat('dd:mm')
    # Scale bar
    ang_sep = (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
    fig.add_colorbar()
    fig.colorbar.set_width(0.15)
    fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, location='top',ticks=[20,20.5,21,21.5,22,22.5,23])
    fig.colorbar.set_font(family='sans_serif',size=text_size)
    fig.add_scalebar(ang_sep.to(u.degree))
    fig.scalebar.set_corner(plot_param['scalebar_pos'])
    fig.scalebar.set(color=label_colour)
    fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
    fig.scalebar.set_font(family='sans_serif',size=text_size)
    fig.tick_labels.set_font(family='sans_serif',size=text_size)
    fig.axis_labels.set_font(family='sans_serif',size=text_size)
    fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                  '{0}\n{1}'.format(region,'log N(H$_2$)'), 
                  relative=True, color=label_colour, 
                  horizontalalignment=plot_param['label_align'],
                  family='sans_serif',size=text_size)
    if include_prot:
        ra_prot, de_prot = get_prot_loc(region)
        if region == 'B18':
            marker_size = 80
            print ra_prot, de_prot
        else: 
            marker_size = 50
        fig.show_markers(ra_prot,de_prot,marker='*',s=marker_size,c='white',edgecolor='black',linewidth=cont_lw*1.5,zorder=4)
        fig.save('figures/{0}_contCol_image_prot.pdf'.format(region),adjust_bbox=True,dpi=100)
    else:
        fig.save( 'figures/{0}_contCol_image_box.pdf'.format(region),adjust_bbox=True)#, bbox_inches='tight')
    fig.close()

def plot_abundance(nh3_cont_fits,nh3_col_hdu,h2_col_hdu,region,plot_pars,maskLim,obsMaskFits):
    text_size = 14
    b18_text_size = 20
    if region == 'B18':
        text_size = b18_text_size
    # Get protostellar locations
    ra_prot, de_prot = get_prot_loc(region)
    # Contour parameters (currently NH3 moment 0)
    cont_color='0.6'
    cont_lw   = 0.6
    cont_levs=2**np.arange( 0,20)*plot_param['w11_step']
    # Calculate abundance
    log_xnh3 = nh3_col_hdu[0].data - np.log10(h2_col_hdu.data)
    log_xnh3_hdu = fits.PrimaryHDU(log_xnh3,nh3_col_hdu[0].header)
    log_xnh3_hdu.writeto('../testing/{0}/parameterMaps/{0}_XNH3_{1}.fits'.format(region,file_extension),clobber=True)
    fig=aplpy.FITSFigure(log_xnh3_hdu,figsize=(plot_param['size_x'], plot_param['size_y']))
    fig.show_colorscale(cmap='YlOrRd_r',vmin=plot_param['xnh3_lim'][0],vmax=plot_param['xnh3_lim'][1])
    #fig.set_nan_color('0.95')
    # Observations mask contour
    fig.show_contour(obsMaskFits,colors='white',levels=1,linewidths=1.5)
    # NH3 moment contours
    # Masking of small (noisy) regions
    selem = np.array([[0,1,0],[1,1,1],[0,1,0]])
    LowestContour = cont_levs[0]*0.5
    w11_hdu = fits.open(nh3_cont_fits)
    map = w11_hdu[0].data
    mask = binary_opening(map > LowestContour, selem)
    MaskedMap = mask*map
    w11_hdu[0].data = MaskedMap
    fig.show_contour(w11_hdu,colors=cont_color,levels=cont_levs,linewidths=cont_lw)
    # Ticks
    fig.ticks.set_color('black')
    fig.tick_labels.set_font(family='sans_serif',size=text_size)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_style('colons')
    fig.tick_labels.set_yformat('dd:mm')
    # Scale bar
    ang_sep = (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
    fig.add_colorbar()
    fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, location='top',
                      ticks=[-10,-9.5,-9,-8.5,-8,-7.5,-7,-6.5])
    fig.colorbar.set_font(family='sans_serif',size=text_size)
    fig.add_scalebar(ang_sep.to(u.degree))
    fig.scalebar.set_font(family='sans_serif',size=text_size)
    fig.scalebar.set_corner(plot_param['scalebar_pos'])
    fig.scalebar.set(color='black')
    fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
    label_colour = 'black'
    fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                  '{0}\n{1}'.format(region,r'$\mathrm{log} \ X(\mathrm{NH}_3)$'), 
                  relative=True, color=label_colour, 
                  horizontalalignment=plot_param['label_align'],
                  family='sans_serif',size=text_size)
    fig.save( 'figures/{0}_xnh3_image.pdf'.format(region),adjust_bbox=True,dpi=200)#, bbox_inches='tight')
    # Add protostars
    fig.show_markers(ra_prot,de_prot,marker='*',s=50,
                     c='white',edgecolors='black',linewidth=0.5,zorder=4)
    fig.save( 'figures/{0}_xnh3_image_prot.pdf'.format(region),adjust_bbox=True,dpi=200)
    fig.close()

def plot_temp_overlay(nh3_temp_hdu,h2_temp_hdu,nh3_cont_fits,region,plot_param,include_prot=True):
    temp_ratio = nh3_temp_hdu[0].data/h2_temp_hdu.data
    temp_ratio[np.where(temp_ratio == 0)] = np.nan
    temp_ratio_hdu = fits.PrimaryHDU(temp_ratio,nh3_temp_hdu[0].header)
    text_size = 14
    b18_text_size = 20    
    # Contour parameters (currently NH3 moment 0)
    cont_color='indigo'
    cont_lw   = 0.6
    cont_levs=2**np.arange( 0,20)*plot_param['w11_step']
    fig=aplpy.FITSFigure(temp_ratio_hdu,figsize=(plot_param['size_x'], plot_param['size_y']))
    fig.show_colorscale(cmap='afmhot',vmin=0.5,vmax=1.5)
    fig.set_nan_color('0.95')
    # Observations mask contour
    fig.show_contour(obsMaskFits,colors='black',levels=1,linewidths=1.5)
    # NH3 moment contours
    # Masking of small (noisy) regions
    selem = np.array([[0,1,0],[1,1,1],[0,1,0]])
    LowestContour = cont_levs[0]*0.5
    w11_hdu = fits.open(nh3_cont_fits)
    map = w11_hdu[0].data
    mask = binary_opening(map > LowestContour, selem)
    MaskedMap = mask*map
    w11_hdu[0].data = MaskedMap
    # Labels
    if region == 'B18':
        text_size = b18_text_size
    label_colour = 'black'
    fig.show_contour(w11_hdu,colors=cont_color,levels=cont_levs,linewidths=cont_lw)
    # Ticks
    fig.ticks.set_color('black')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_style('colons')
    fig.tick_labels.set_yformat('dd:mm')
    # Scale bar
    ang_sep = (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
    fig.add_colorbar()
    fig.colorbar.set_width(0.15)
    fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, location='top')
    fig.colorbar.set_font(family='sans_serif',size=text_size)
    fig.add_scalebar(ang_sep.to(u.degree))
    fig.scalebar.set_corner(plot_param['scalebar_pos'])
    fig.scalebar.set(color='black')
    fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
    fig.scalebar.set_font(family='sans_serif',size=text_size)
    fig.tick_labels.set_font(family='sans_serif',size=text_size)
    fig.axis_labels.set_font(family='sans_serif',size=text_size)
    fig.add_label(plot_param['label_xpos'], plot_param['label_ypos'], 
                  '{0}\n{1}'.format(region,'$T_K \ / \ T_d$'), 
                  relative=True, color=label_colour, 
                  horizontalalignment=plot_param['label_align'],
                  family='sans_serif',size=text_size)
    if include_prot:
        ra_prot, de_prot = get_prot_loc(region)
        if region == 'B18':
            marker_size = 70
        else: 
            marker_size = 50
        fig.show_markers(ra_prot,de_prot,marker='*',s=marker_size,c='white',edgecolor='black',linewidth=cont_lw*1.5,zorder=4)
        fig.save('figures/{0}_tempRatio_image_prot.pdf'.format(region),adjust_bbox=True,dpi=120)
    else:
        fig.save( 'figures/{0}_tempRatio_image.pdf'.format(region),adjust_bbox=True,dpi=120)#, bbox_inches='tight')
    fig.close()    

herDir  = '../otherData/herschel_ayushi/'

region_list  = ['L1688','NGC1333','B18','OrionA']
herFile_list = ['OphL1688','perseus','Tau_B18','orionA-N']
eTkin_lim = 1.
eTex_lim = 2.
file_extension = 'DR1_rebase3'
maskLimArray = [0.25,0.5,0.5,0.2]

for region_i in range(len(region_list)):
    region = region_list[region_i]
    plot_param = plottingDictionary[region]
    maskLim = maskLimArray[region_i]
    her_file_name = herFile_list[region_i]
    nh3ImFits  = '../testing/{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension)
    nh3RmsFits = '../testing/{0}/{0}_NH3_11_{1}_rms_QA_trim.fits'.format(region,file_extension)
    nh3ColFits = '../testing/{0}/parameterMaps/{0}_N_NH3_{1}_flag.fits'.format(region,file_extension)
    nh3TempFits = '../testing/{0}/parameterMaps/{0}_Tkin_{1}_flag.fits'.format(region,file_extension)
    nh3eTempFits = '../testing/{0}/parameterMaps/{0}_eTkin_{1}_flag.fits'.format(region,file_extension)
    nh3eTexFits = '../testing/{0}/parameterMaps/{0}_eTex_{1}_flag.fits'.format(region,file_extension)
    nh3PropMaps = '../testing/{0}/{0}_parameter_maps_{1}_trim.fits'.format(region,file_extension)
    obsMaskFits = '../testing/{0}/{0}_NH3_11_{1}_obsMask.fits'.format(region,file_extension)
    herImFits  = herDir+'{0}/Herschel_data/{1}-500.offset.fits'.format(region,her_file_name)
    herColFits = herDir+'{0}/Colden_temp/{1}_colden_masked.fits'.format(region,her_file_name)
    herTdFits  = herDir+'{0}/Colden_temp/{1}_temp_masked.fits'.format(region,her_file_name)
    # Compare line brightness, continuum flux density
    nh3_fits = fits.open(nh3ImFits)
    h2_regrid_fits = regrid_h2(nh3ImFits,herImFits)
    # Compare H2, NH3 column densities
    nh3_col = fits.open(nh3ColFits)
    nh2_regrid_fits = regrid_h2(nh3ColFits,herColFits)
    # Write out regridded fits file if doesn't exist
    if not os.path.isfile('{0}_NH2_regrid.fits'.format(region)):
        new_hdulist = fits.HDUList([nh2_regrid_fits])
        new_hdulist.writeto('{0}_NH2_regrid.fits'.format(region))
    # Mask where no good fits (data == 0)
    nh3_col[0].data[(nh3_col[0].data == 0)] = np.nan
    # Mask where uncertainties on Tex are high
    nh3eTex = fits.getdata(nh3eTexFits)
    nh3_col[0].data[(nh3eTex > eTex_lim)] = np.nan
    # Plot NH3 m0 contours over H2 column, X(NH3)
    nh2_log_fits = log10_h2(nh2_regrid_fits)
    plot_column_overlay(nh3ImFits,nh2_log_fits,region,plot_param,maskLim,obsMaskFits,include_prot=True)
    plot_abundance(nh3ImFits,nh3_col,nh2_regrid_fits,region,plot_param,maskLim,obsMaskFits)
    # Compare H2, NH3 temperatures
    # Add mask based on uncertainty in Tk
    # Mask where no good fits (Tkin == 0)
    nh3_temp = fits.open(nh3TempFits)
    nh3_etemp = fits.open(nh3eTempFits)
    # Plot images of Tk/Td ratio
    h2_temp_regrid_fits = regrid_h2(nh3TempFits,herTdFits)
    plot_temp_overlay(nh3_temp,h2_temp_regrid_fits,nh3ImFits,region,plot_param,include_prot=True)
    # Plot images of H2 with NH3 contours
    vmin = np.nanmin(h2_regrid_fits.data)
    vmax = np.nanmax(h2_regrid_fits.data)
    plot_dust_overlay(nh3ImFits,h2_regrid_fits,region,plot_param,vmin,vmax,maskLim,obsMaskFits)
