import astropy.units as u
import numpy as np

L1688 = {
    'vmin' : -25,
    'vmax' : 40,
    #'max_amp' : 20000,
    #'y_offsets' : np.array([7000,3500,2400,1600,800,0]),
    'min_amp' : -0.03,
    'max_amp' : 0.9,
    'y_offsets' : np.array([0.5,0.4,0.3,0.2,0.1,0]),
    'label_offsets': np.array([0.1,0.03,0.02,0.02,0.02,0.02]),
    'mask33_chans' : np.array([160,216]) #6,10
    }

B18 = {
    'vmin' : -25,
    'vmax' : 40,
#    'min_amp' : -0.25,
#    'max_amp' : 10000,
#    'y_offsets' : np.array([5000,4200,3500,1300,500,0]),
    'min_amp' : -0.03,
    'max_amp' : 0.4,
    'y_offsets' : np.array([0.3,0.25,0.2,0.1,0.05,0]),
    'label_offsets': np.array([0.03,0.02,0.02,0.01,0.01,0.01]),
#    'label_offsets': np.array([1400,400,250,100,100,100]),
    'mask33_chans' : np.array([109,207]) # 18,25
    }

NGC1333 = {
    'vmin' : -25,
    'vmax' : 40,
#    'min_amp' : -0.25,
#    'max_amp' : 8000,
#    'y_offsets' : np.array([3500,2500,1800,1000,500,0]),
#    'label_offsets': np.array([1400,300,250,100,100,100]),
    'min_amp' : -0.03,
    'max_amp' : 0.45,
    'y_offsets' : np.array([0.25,0.2,0.15,0.1,0.05,0]),
    'label_offsets': np.array([0.06,0.02,0.02,0.02,0.02,0.02]),
    'mask33_chans' : np.array([346,396]) #24,27
    }

OrionA = {
    'vmin' : -25,
    'vmax' : 40,
#    'min_amp' : -0.25,
#    'max_amp' : 10500,
#    'y_offsets' : np.array([5500,3500,2400,1600,800,0]),
#    'label_offsets': np.array([1600,500,400,100,100,100]),
    'min_amp' : -0.03,
    'max_amp' : 0.45,
    'y_offsets' : np.array([0.3,0.225,0.15,0.1,0.05,0]),
    'label_offsets': np.array([0.1,0.03,0.02,0.01,0.01,0.01]),
    'mask33_chans' : np.array([130,192]) #14.6,19
    }

plottingDictionary = {"L1688" : L1688, "B18" : B18, "NGC1333" : NGC1333, "OrionA" : OrionA}
