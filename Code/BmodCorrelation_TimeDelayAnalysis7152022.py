# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 11:07:01 2022

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

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y


import smooth as sm
def iter_smooth(array,loops=6,window_len=3):
    for l in np.arange(loops):
        array = sm.smooth(array,window_len=window_len)
    return array

# Sample rate and desired cutoff frequencies (in Hz).
fs = 125e6
lowcut = 50e3
highcut = 2e6

directory='C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'
datafilename='7152022_3p5kV_51shots.h5'
data=load_hdf5(directory+datafilename,verbose=True)


time_s = data['time']['time_s'][:]
timeB_s = data['time']['timeB_s'][:]
time_us = data['time']['time_us'][:]
timeB_us = data['time']['timeB_us'][:]
analysis_start_time = 60
analysis_end_time = 160
start_time_index = iff.tindex_min(analysis_start_time,timeB_us[0,:])
end_time_index = iff.tindex_min(analysis_end_time,timeB_us[0,:])
timerange_limit = 3e-6#s
port_sep = 0.0254#m

numshots=94
direction_list=['r','theta','z']
num_pairs = 3
directions = 3
#tde_pairs=[['pos5','pos7'],['pos19','pos21'],['pos33','pos35']]
tde_pairs=[['pos19','pos21']]
tde_pairs=[['pos5','pos7']]

delaytimes = np.zeros([num_pairs, numshots])
delayindex = np.zeros([num_pairs, numshots])

plotshots=1
plt.rc('axes',linewidth=0.75)
plt.rc('xtick.major',width=0.75)
plt.rc('ytick.major',width=0.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)
plt.rc('lines',markersize=2.5,markeredgewidth=0.0)

fig=plt.figure(num=1,figsize=(6,3),dpi=300,facecolor='w',edgecolor='k')

left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.90      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.05   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
#ax=plt.axes([left,bottom,right-left,top-bottom])
ax=plt.subplot(1,1,1)

# Function to print mouse click event coordinates
def onclick(event):
    if event.key=='shift':
        print(event.xdata)

for shot in np.arange(46,52):
    tde_pair_index=0
    for tde_pair in tde_pairs:
        
        data1r=data[tde_pair[0]]['Bfilt']['r'][shot,:]
        data1t=data[tde_pair[0]]['Bfilt']['theta'][shot,:]
        data1z=data[tde_pair[0]]['Bfilt']['z'][shot,:]
        data1mod=np.sqrt(data1r**2+data1t**2+data1z**2)
        data1mod_max=np.max(data1mod[start_time_index:end_time_index])
        data1mod_norm=data1mod/data1mod_max
        data1mod_sm=iter_smooth(data1mod_norm,loops=3,window_len=11)
        
        data2r=data[tde_pair[1]]['Bfilt']['r'][shot,:]
        data2t=data[tde_pair[1]]['Bfilt']['theta'][shot,:]
        data2z=data[tde_pair[1]]['Bfilt']['z'][shot,:]
        data2mod=np.sqrt(data2r**2+data2t**2+data2z**2)
        data2mod_max=np.max(data2mod[start_time_index:end_time_index])
        data2mod_norm=data2mod/data2mod_max
        data2mod_sm=iter_smooth(data2mod_norm,loops=3,window_len=11)

        #data1mod_filt=butter_bandpass_filter(data1mod_norm,lowcut,highcut,fs,order=9)
        #data2mod_filt=butter_bandpass_filter(data2mod_norm,lowcut,highcut,fs,order=9)
        
        #correlation
        t=timeB_us[0, start_time_index:end_time_index]
        d1=data1mod_norm[start_time_index:end_time_index]
        d2=data2mod_norm[start_time_index:end_time_index]
        tau,corr=gc.get_corr(t,d2,d1,normalized=False)
        
        
        #use filtered
        #d1=data1_filt[start_time_index:end_time_index]
        #d2=data2_filt[start_time_index:end_time_index]
        if plotshots:
            
            plt.plot(np.array(timeB_us[0,:]),data1mod_sm, linewidth=0.5)
            plt.plot(np.array(timeB_us[0,:]),data2mod_sm, linewidth=0.5)
            plt.title(tde_pair[0]+tde_pair[1]+' Shot '+str(shot) + 'B_mod')
            plt.xlim(0.0,40.0)
            plt.ylim(0,0.5)
            plt.xlabel(r'Time ($\mu$s)')
            plt.ylabel('B Modulus')
            """
            # Bind the button_press_event with the onclick() method
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()    
            plt.ginput(3,timeout=0)
            plt.clf()
            
            plt.plot(tau,corr, linewidth=0.5)
            plt.title(tde_pair[0]+tde_pair[1]+' Shot '+str(shot) + 'Correlation')
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()    
            plt.ginput(6,timeout=0)
            plt.clf()"""
            
           #savefilename=tde_pair[0]+tde_pair[1]+'_shot_'+str(shot+1).zfill(2)+'mod.png'
            #savefile = savedirectory+savefilename
            #plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
            #plt.clf()
        #use unfiltered
        #d1=data1_norm[start_time_index:end_time_index]
        #d2=data2_norm[start_time_index:end_time_index]

        #t=timeB_s[start_time_index:end_time_index]
        #dt=timeB_s[1]-timeB_s[0]
        #tau,corr=gc.get_corr(t,data2mod_norm,data1mod_norm,normalized=False)
        
        #index_at_zero_time=iff.tindex_min(0,tau)
        #indexsize_of_timerange_limit = int(np.round(timerange_limit/dt))

        #find time (in seconds) of max in correlation function within timerange limit
        #delayindex[tde_pair_index,shot]=np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))
        #delay=tau[index_at_zero_time-indexsize_of_timerange_limit+np.argmax(np.abs(corr[index_at_zero_time-indexsize_of_timerange_limit:index_at_zero_time+indexsize_of_timerange_limit]))]
        #delaytimes[tde_pair_index,shot]=delay
    tde_pair_index+=1
            
