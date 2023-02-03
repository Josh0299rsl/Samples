# -*- coding: utf-8 -*-
"""
Created on Mon May 16 16:52:22 2022

@author: Josh0
"""

import scipy.io as spio
import matplotlib.pylab as plt
import numpy as np
import os
from load_hdf5 import load_hdf5
import spectrum_wwind as spec
import indexfinderfuncs as iff
import get_corr as gc
import scipy.integrate as sp
from scipy.interpolate import interp1d

from scipy.signal import butter, sosfiltfilt, sosfreqz
from scipy.signal import savgol_filter as sv
from scipy.stats import norm

fs = 125e6
lowcut = 50e3
highcut = 2e6

plt.style.use('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\__pycache__\\mplstyle_tallplots.py')
directory='C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'
datafilename='3172022_253shots_checkpy_forbadshots.h5'
data=load_hdf5(directory+datafilename,verbose=True)


time_s = data['time']['time_s'][:]
timeB_s = data['time']['timeB_s'][:]
time_us = data['time']['time_us'][:]
timeB_us = data['time']['timeB_us'][:]
analysis_start_time = 50
analysis_end_time = 150
start_time_index = iff.tindex_min(analysis_start_time,timeB_us[0,:])
end_time_index = iff.tindex_min(analysis_end_time,timeB_us[0,:])
timerange_limit = 3e-6#s
port_sep = 0.0254#m

numshots=81
direction_list=['r','theta','z']
num_pairs = 2
directions = 3
tde_pairs = [['pos19','pos21']]

#BmagMean_inner = np.zeros()
#Instead of having arrays with just numbers of shots, I have to do full size and remove zeros.
BmagMean_fullnoz19 = np.zeros(253) #(122) Same for all below
BmagMean_halfnoz19 = np.zeros(253) #(69) Same for all below

BmagMean_fullnoz21 = np.zeros(253)
BmagMean_halfnoz21 = np.zeros(253)

valfven_full19 = np.zeros(253)
valfven_half19 = np.zeros(253)

valfven_full21 = np.zeros(253)
valfven_half21 = np.zeros(253)

alfven_const = 0.004587

for shot in np.arange(57, 180):
   
    datar19=data['pos19']['B']['r'][shot,start_time_index:end_time_index]
    datat19=data['pos19']['B']['theta'][shot,start_time_index:end_time_index]
    dataz19=data['pos19']['B']['z'][shot,start_time_index:end_time_index]
    datamod19=np.sqrt(datar19**2+datat19**2+dataz19**2)
    
    datar21=data['pos21']['B']['r'][shot,start_time_index:end_time_index]
    datat21=data['pos21']['B']['theta'][shot,start_time_index:end_time_index]
    dataz21=data['pos21']['B']['z'][shot,start_time_index:end_time_index]
    datamod21=np.sqrt(datar21**2+datat21**2+dataz21**2)

    BmagMean_fullnoz19[shot] = np.mean(datamod19[shot])
    BmagMean_fullnoz21[shot] = np.mean(datamod21[shot])
    valfven_full19[shot] = BmagMean_fullnoz19[shot]/alfven_const*1e-3
    valfven_full21[shot] = BmagMean_fullnoz21[shot]/alfven_const*1e-3

for shot in np.arange(184, 253):
    
    datar19=data['pos19']['B']['r'][shot,start_time_index:end_time_index]
    datat19=data['pos19']['B']['theta'][shot,start_time_index:end_time_index]
    dataz19=data['pos19']['B']['z'][shot,start_time_index:end_time_index]
    datamod19_half=np.sqrt(datar19**2+datat19**2+dataz19**2)
    
    datar21=data['pos21']['B']['r'][shot,start_time_index:end_time_index]
    datat21=data['pos21']['B']['theta'][shot,start_time_index:end_time_index]
    dataz21=data['pos21']['B']['z'][shot,start_time_index:end_time_index]
    datamod21_half=np.sqrt(datar21**2+datat21**2+dataz21**2)

    BmagMean_halfnoz19[shot] = np.mean(datamod19_half[shot])
    BmagMean_halfnoz21[shot] = np.mean(datamod21_half[shot])
    valfven_half19[shot] = BmagMean_halfnoz19[shot]/alfven_const*1e-3
    valfven_half21[shot] = BmagMean_halfnoz21[shot]/alfven_const*1e-3


bins = np.array(range(50, 1000, 30))

cut_valfven_full19 = valfven_full19[(valfven_full19 >= 100) & (valfven_full19 <= 1000)]
cut_valfven_full21 = valfven_full21[(valfven_full21 >= 100) & (valfven_full21 <= 1000)]
cut_valfven_half19 = valfven_half19[(valfven_half19 >= 100) & (valfven_half19 <= 1000)]
cut_valfven_half21 = valfven_half21[(valfven_half21 >= 100) & (valfven_half21 <= 1000)]

(mu19_full, sigma19_full) = norm.fit(cut_valfven_full19)
(mu21_full, sigma21_full) = norm.fit(cut_valfven_full21)

(mu19_half, sigma19_half) = norm.fit(cut_valfven_half19)
(mu21_half, sigma21_half) = norm.fit(cut_valfven_half21)

mu19_full = round(mu19_full)
mu21_full = round(mu21_full)
mu19_half = round(mu19_half)
mu21_half = round(mu21_half)

fig, (ax1, ax2) = plt.subplots(2)
ax1.hist(valfven_half19, bins=bins, color='red', edgecolor='black', label=r'Pos19, Half Noz, $\mu_{19}$='+str(mu19_half)+' km/s', align='left')
ax2.hist(valfven_half21, bins=bins, color='teal', edgecolor='black', label=r'Pos21, Half Noz, $\mu_{21}$='+str(mu21_half)+' km/s', align='left')
ax1.set_xlabel(r'V$_a$ (km/s)')
ax2.set_xlabel(r'V$_a$ (km/s)')
ax1.set_ylabel('Counts')
ax2.set_ylabel('Counts')
ax1.legend(loc='best', frameon=False)
ax2.legend(loc='best', frameon=False)
#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\2232022_3172022\\MagStructures_BulkVelocity\\AlfvenSpeedStudies\\Pos1921HalfNozzle_alfvenHistogram.png', dpi=300)













