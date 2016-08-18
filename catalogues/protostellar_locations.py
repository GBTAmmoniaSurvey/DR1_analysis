import numpy as np
'''
Snippet of code to parse the protostellar catalogue files in this directory. 
Catalogues have been edited to more easily be read in (_edit.txt)
'''
def get_prot_loc(region_name):
    # C2D doesn't include Orion, Taurus
    c2d_data_file = 'c2d_dunham_table2_edit.txt'
    orion_file    = 'megeath_orion_table4_edit.txt'
    taurus1_file  = 'rebull_taurus_table6_edit.txt'
    taurus2_file  = 'rebull_taurus_table7_edit.txt'
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
        ra_deg = 15. * (protData['ra1']+(protData['ra2'] + protData['ra3']/60.)/60.)
        de_deg = np.abs(protData['de1'])+(protData['de2'] + protData['de3']/60.)/60.
        de_deg[protData['de1'] < 0] = de_deg[protData['de1']<0] * (-1.)
        ra_deg_class0I_2 = ra_deg[np.logical_or(protData['class'] == 'I',
                                                protData['class'] == 'flat')]
        de_deg_class0I_2 = de_deg[np.logical_or(protData['class'] == 'I',
                                                protData['class'] == 'flat')]
        ra_deg_class0I = np.append(ra_deg_class0I_1,ra_deg_class0I_2)
        de_deg_class0I = np.append(de_deg_class0I_1,de_deg_class0I_2)
    
    return ra_deg_class0I, de_deg_class0I
