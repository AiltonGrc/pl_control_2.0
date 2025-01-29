# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 17:07:37 2022

@author: Gabriel Undeutsch
"""

import spe_loader as sl
import matplotlib.pyplot as plt


def get_data_from_spe(file):
    # call the spe-file using the spe_loader library
    spe_files = sl.load_from_file(file)
    
    # get wavelength and intensity out of the spe-file
    wavelength = spe_files.wavelength
    intensity = []
    
    for num_map in range(spe_files.nframes):
        intensity.append(spe_files.data[num_map][0][0])
    
    return wavelength, intensity

# file = r'I:/public/NANOSCALE SEMICONDUCTOR GROUP/1. DATA/BIG-LAB/2022/07/11/SA323_QD4_green_1s#Pol(deg)_0_360.spe'
# file = r'I:/public/NANOSCALE SEMICONDUCTOR GROUP/1. DATA/SMALL-LAB/2022/March/30/SA0484_QD1_50umslit_1800grating_greenlaser_FEL600_1s_7K#Pol(deg)_0_360.spe'
# l,i = get_data_from_spe(file)
#
# plt.plot(l,i[22])


