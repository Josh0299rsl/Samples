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
            
            #savedirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\timeseries\\probepairs\\'+tde_pair[0]+tde_pair[1]+'\\mod\\'
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
            
#velocities = (port_sep/delaytimes)/1000.0#km/s
#mean_velocities = np.mean(velocities,axis=1)

#np.savez(directory+'mean_vels_01122022_60t160_50to500kHzfilt',mean_velocities=mean_velocities)
"""
plt.figure(57)
plt.plot(velocities[0,0,:],'o')
plt.plot(velocities[0,1,:],'o')
plt.plot(velocities[0,2,:],'o')
plt.plot(mean_velocities[0,:],'x',color='red')
plt.figure(1921)
plt.plot(velocities[1,0,:],'o')
plt.plot(velocities[1,1,:],'o')
plt.plot(velocities[1,2,:],'o')
plt.plot(mean_velocities[1,:],'x',color='red')
plt.figure(3335)
plt.plot(velocities[2,0,:],'o')
plt.plot(velocities[2,1,:],'o')
plt.plot(velocities[2,2,:],'o')
plt.plot(mean_velocities[2,:],'x',color='red')
#########plot velocity distribution ###############
plt.rc('axes',linewidth=2.0)
plt.rc('xtick.major',width=2.0)
plt.rc('ytick.major',width=2.0)
plt.rc('xtick.minor',width=2.0)
plt.rc('ytick.minor',width=2.0)
plt.rc('lines',markersize=8,markeredgewidth=0.0,linewidth=1.0)
fig=plt.figure(num=571,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(mean_velocities[0,:], bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
#plt.xlim(0,198)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
#plt.ylim(0,5)
#plt.text(0.50,0.92,r'Mean: '+mean_vel_str+'$\pm$'+std_vel_str+'km/s',transform=ax1.transAxes,fontsize=12)
fig=plt.figure(num=19211,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(mean_velocities[1,:], bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
#plt.xlim(0,198)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
#plt.ylim(0,5)
#plt.text(0.50,0.92,r'Mean: '+mean_vel_str+'$\pm$'+std_vel_str+'km/s',transform=ax1.transAxes,fontsize=12)
fig=plt.figure(num=33351,figsize=(5,4),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.2  # the bottom of the subplots of the figure
top = 0.96      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.25   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax1=plt.subplot(1,1,1)
plt.hist(mean_velocities[2,:], bins=20,range=(0,100))  # arguments are passed to np.histogram
plt.xticks(fontsize=12)
plt.xlabel(r'Bulk Vel. [km/s]',fontsize=16)
#plt.xlim(0,198)
#plt.yticks(np.array([0,2,4,6,8,10]),[0,2,4,6,8,10],fontsize=12)
plt.ylabel('Count',fontsize=16)
#plt.xlim(50,82)
#plt.ylim(0,5)
#plt.text(0.50,0.92,r'Mean: '+mean_vel_str+'$\pm$'+std_vel_str+'km/s',transform=ax1.transAxes,fontsize=12)
#b5r=data['mag_probe']['positions']['probe5']['r']['b'][shot,:]
#b5r_max=np.max(np.abs(b5r))
#b5r_norm=b5r/b5r_max
#b7r=data['mag_probe']['positions']['probe7']['r']['b'][shot,:]
#b7r_max=np.max(np.abs(b7r))
#b7r_norm=b7r/b7r_max
#b5t=data['mag_probe']['positions']['probe5']['t']['b'][shot,:]
#b5t_max=np.max(np.abs(b5t))
#b5t_norm=b5t/b5t_max
#b7t=data['mag_probe']['positions']['probe7']['t']['b'][shot,:]
#b7t_max=np.max(np.abs(b7t))
#b7t_norm=b7t/b7t_max
#b5z=data['mag_probe']['positions']['probe5']['z']['b'][shot,:]
#b5z_max=np.max(np.abs(b5z))
#b5z_norm=b5z/b5z_max
#b7z=data['mag_probe']['positions']['probe7']['z']['b'][shot,:]
#b7z_max=np.max(np.abs(b7z))
#b7z_norm=b7z/b7z_max
#b5rfilt=butter_bandpass_filter(b5r_norm,lowcut,highcut,fs,order=9)
#b7rfilt=butter_bandpass_filter(b7r_norm,lowcut,highcut,fs,order=9)
#b5tfilt=butter_bandpass_filter(b5t_norm,lowcut,highcut,fs,order=9)
#b7tfilt=butter_bandpass_filter(b7t_norm,lowcut,highcut,fs,order=9)
#b5zfilt=butter_bandpass_filter(b5z_norm,lowcut,highcut,fs,order=9)
#b7zfilt=butter_bandpass_filter(b7z_norm,lowcut,highcut,fs,order=9)
#plt.figure(1)
#plt.clf()
#plt.plot(timeB_us,b5r)
#plt.plot(timeB_us,b5rfilt)
#plt.plot(timeB_us,b7rfilt)
# Plot the frequency response for a few different orders.
#plt.figure(1)
#plt.clf()
#for order in [3, 6, 9]:
#    sos = butter_bandpass(lowcut, highcut, fs, order=order)
#    w, h = sosfreqz(sos, worN=2000)
#    plt.semilogx((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
#d5r=b5rfilt[start_time_index:end_time_index]
#d7r=b7rfilt[start_time_index:end_time_index]
#d5t=b5tfilt[start_time_index:end_time_index]
#d7t=b7tfilt[start_time_index:end_time_index]
#d5z=b5zfilt[start_time_index:end_time_index]
#d7z=b7zfilt[start_time_index:end_time_index]
#t=timeB_us[start_time_index:end_time_index]
#tau57r,corr57r=gc.get_corr(t,d7r,d5r,normalized=False)
#tau57t,corr57t=gc.get_corr(t,d7t,d5t,normalized=False)
#tau57z,corr57z=gc.get_corr(t,d7z,d5z,normalized=False)
#plt.plot(tau57r,corr57r)
#plt.plot(tau57t,corr57t)
#plt.plot(tau57z,corr57z)
"""

"""
plt.rc('axes',linewidth=0.75)
plt.rc('xtick.major',width=0.75)
plt.rc('ytick.major',width=0.75)
plt.rc('xtick.minor',width=0.75)
plt.rc('ytick.minor',width=0.75)
plt.rc('lines',markersize=2.5,markeredgewidth=0.0)
fig=plt.figure(num=1,figsize=(6,3),dpi=300,facecolor='w',edgecolor='k')
left  = 0.15  # the left side of the subplots of the figure
right = 0.97    # the right side of the subplots of the figure
bottom = 0.05   # the bottom of the subplots of the figure
top = 0.97      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.05   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax=plt.axes([left,bottom,right-left,top-bottom])
for shot in np.arange(111,112):
    pico1=spio.loadmat(directory+'Pico1\\20220112-0001 ('+str(shot)+').mat')
    pico2=spio.loadmat(directory+'Pico2\\20220112-0001 ('+str(shot)+').mat')
    pico3=spio.loadmat(directory+'Pico3\\20220112-0001 ('+str(shot)+').mat')
    pico4=spio.loadmat(directory+'Pico4\\20220112-0001 ('+str(shot)+').mat')
    pico5=spio.loadmat(directory+'Pico5\\20220112-0001 ('+str(shot)+').mat')
    
    ax=plt.subplot(3,2,1)
    plt.plot(pico1['A'])
    plt.plot(pico1['B'])
    plt.plot(pico1['C'])
    plt.plot(pico1['D'])
    
    ax=plt.subplot(3,2,2)
    plt.plot(pico2['A'])
    plt.plot(pico2['B'])
    plt.plot(pico2['C'])
    plt.plot(pico2['D'])
    
    ax=plt.subplot(3,2,3)
    plt.plot(pico3['A'])
    plt.plot(pico3['B'])
    plt.plot(pico3['C'])
    plt.plot(pico3['D'])
    
    ax=plt.subplot(3,2,4)
    plt.plot(pico4['A'])
    plt.plot(pico4['B'])
    plt.plot(pico4['C'])
    plt.plot(pico4['D'])
    
    ax=plt.subplot(3,2,5)
    plt.plot(pico5['A'])
    plt.plot(pico5['B'])
    plt.plot(pico5['C'])
    plt.plot(pico5['D'])
    
    save_dir = 'C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\QuickPlots\\'
    filename = 'Shot'+str(shot)+'.png'
    savefile = os.path.normpath(save_dir+filename)
    plt.savefig(savefile,dpi=300,facecolor='w',edgecolor='k')
    plt.clf()
   """ 