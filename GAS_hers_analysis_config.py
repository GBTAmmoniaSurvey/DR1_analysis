import astropy.units as u
import numpy as np

L1688={'size_x' : 9, 'size_y' : 7, 
    'scalebar_size' : 0.1*u.pc, 'scalebar_pos' : 'bottom right',
    'distance' : 137.3*u.pc, 'beam_pos' : 'bottom left',
    'label_xpos' : 0.025, 'label_ypos' : 0.1,
    'label_align' : 'left',
    'rms_min' : np.array([0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05]),
    'rms_max' : np.array([0.75, 0.75, 0.75, 1.0, 0.75, 0.75, 0.75]),
    'mom0_min' : np.array([-0.14, -0.12, -0.15, -0.35, -0.14, -0.14, -0.18]), 
    'mom0_max' : np.array([11.0,   1.8,   0.4,  1.4,   0.70,  0.70,  0.90]),
    'mom0_stretch' : 'log',
    'pmin_list': np.array([-5,0,5,2.73,13]),
    'pmax_list': np.array([20,2,30,15,15.5]),
    'p_step': np.array([0.4,0.2,4,2,0.5]),
    'w11_step' : 0.3*3,
    # Herschel comparison parameters
    'vmid_scale':2,
    'nh2_grid_size':70,
    'temp_grid_size':60,
    'flux_grid_size':200,
    'her_flux_max':1000.,
    'nh3_flux_max':15.,
    'nh3_col_lim':np.array([13.2,14.7]),
    'h2_col_lim':np.array([21.2,23.1]),
    'cumulative_bins':20,
    'temp_lim':np.array([5.,30.]),
    'xnh3_lim':np.array([-9,-7])}

B18={'size_x' : 19, 'size_y' : 7, 
    'scalebar_size' : 0.1*u.pc, 'scalebar_pos' : 'bottom right',
    'distance' : 135.*u.pc, 'beam_pos' : 'top left',
    'label_xpos' : 0.025, 'label_ypos' : 0.65,
    'label_align' : 'left',
    'rms_min' : np.array([0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05]),
    'rms_max' : np.array([0.75, 0.75, 0.75, 1.0, 0.75, 0.75, 0.75]),
    'mom0_min' : np.array([-0.17,-0.1,-0.1,-0.3,-0.1,-0.2,-0.2]),
    'mom0_max' : np.array([7.5, 0.5, 0.5, 1.8, 1.1, 1.0, 1.5]),
    'mom0_stretch' : 'linear',
    'pmin_list': np.array([4.5,0,5,2.73,13]),
    'pmax_list': np.array([7.5,0.5,15,15,14.5]),
    'p_step': np.array([1,0.1,2,1,0.5]),
    'w11_step' : 0.3*3,
    # Herschel comparison parameters
    'vmid_scale':10,
    'nh2_grid_size':45,
    'temp_grid_size':50,
    'flux_grid_size':200,
    'her_flux_max':210.,
    'nh3_flux_max':10.,
    'nh3_col_lim':np.array([13.2,14.7]),
    'h2_col_lim':np.array([21.3,22.5]),
    'cumulative_bins':20,
    'temp_lim':np.array([5.,20.]),
    'xnh3_lim':np.array([-9,-7])}

NGC1333={'size_x' : 8, 'size_y' : 13.3, 
    'scalebar_size' : 0.1*u.pc, 'scalebar_pos' : 'top right',
    'distance' : 260.*u.pc, 'beam_pos' : 'bottom left',
    'label_xpos' : 0.025, 'label_ypos' : 0.075,
    'label_align' : 'left',
    'rms_min' : np.array([0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05]),
    'rms_max' : np.array([0.75, 0.75, 0.75, 1.0, 0.75, 0.75, 0.75]),
    'mom0_min' : np.array([-0.1, -0.2, -0.1, -0.2, -0.2, -0.2, -0.2]),
    'mom0_max' : np.array([ 15.0,  2.0 ,  0.5,  1.8,  0.5,  0.8,  1.5]),
    'mom0_stretch' : 'log',
    'pmin_list': np.array([-5,0,5,2.73,13]),
    'pmax_list': np.array([20,1,25,15,15]),
    'p_step': np.array([1,0.2,4,1,0.5]),
    'w11_step' : 0.3*2.,
    # Herschel comparison parameters
    'vmid_scale':2,
    'nh2_grid_size':75,
    'temp_grid_size':60,
    'flux_grid_size':200,
    'her_flux_max':800.,
    'nh3_flux_max':15.,
    'nh3_col_lim':np.array([13.2,15.]),
    'h2_col_lim':np.array([21.2,23]),
    'cumulative_bins':20,
    'temp_lim':np.array([5.,30.]),
    'xnh3_lim':np.array([-9,-7])}

OrionA={'size_x' : 4.5, 'size_y' : 11.25, 
    'scalebar_size' : 0.1*u.pc, 'scalebar_pos' : 'bottom right',
    'distance' : 414.*u.pc, 'beam_pos' : 'top left',
    'label_xpos' : 0.025, 'label_ypos' : 0.925,
    'label_align' : 'left',
    'rms_min' : np.array([0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05]),
    'rms_max' : np.array([0.75, 0.75, 0.75, 1.0, 0.75, 0.75, 0.75]),
    'mom0_min' : np.array([-0.1, -0.2, -0.2, -0.1, -0.1, -0.1, -0.2]),
    'mom0_max' : np.array([ 45.0,  5.0 ,  3.0,  1.,  0.5,  0.8,  1.]),
    'mom0_stretch' : 'log',
    'pmin_list': np.array([-5,0,5,2.73,13]),
    'pmax_list': np.array([13,2,40,15,15.5]),
    'p_step': np.array([1,0.4,5,3,0.5]),
    'w11_step' : 0.3*3,
    # Herschel comparison parameters
    'vmid_scale':2,
    'nh2_grid_size':100,
    'temp_grid_size':80,
    'flux_grid_size':200,
    'her_flux_max':12000.,
    'nh3_flux_max':20.,
    'nh3_col_lim':np.array([13.2,15.6]),
    'h2_col_lim':np.array([21.2,24.]),
    'cumulative_bins':20,
    'temp_lim':np.array([5.,50.]),
    'xnh3_lim':np.array([-9.5,-6.5])}

plottingDictionary = {"L1688" : L1688, "B18" : B18, "NGC1333" : NGC1333, "OrionA" : OrionA}

