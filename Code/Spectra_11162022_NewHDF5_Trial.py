# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:42:03 2022

@author: Josh0
"""


import matplotlib.pylab as plt
from load_hdf5 import load_hdf5
import numpy as np
import spectrum_wwind as spec
import indexfinderfuncs as iff

plt.style.use('C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Code\\mplstyle_presentationplots.py')

# Load the file
data_directory_location='C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'
data_directory_location2 = 'C:\\Users\\Josh0\\OneDrive - brynmawr.edu\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Processed Data Files\\2022\\'

#datafilename = '03102020_25shots.h5'
datafilename = '7202022_3p5kV_203shots_Mixof225V_250V_or5ms_Nozzle.h5' #Fast!
datafilename2 = '9232022_3p5kV_NozzleSweep_144shots.h5' #Slow!

data = load_hdf5(data_directory_location + datafilename, verbose=True)
data2 = load_hdf5(data_directory_location2 + datafilename2, verbose=True)

timeB_us = data['time']['timeB_us'][:]
timeB_s = timeB_us*1e-6
time_s = data['time']['time_s'][:]

timeB_us2 = data2['time']['timeB_us'][:]
timeB_s2 = timeB_us2*1e-6
time_s2 = data2['time']['time_s'][:]

poslist = ['pos5', 'pos7', 'pos19','pos21','pos33','pos35']
poslist2 = ['pos5','pos6','pos7','pos8']


start_time = 50
end_time = 150
time_range = [start_time,end_time]

start_time1 = 30
end_time1 = 60
time_range1 = [start_time1,end_time1]
start1_index = iff.tindex_min(start_time1, timeB_us)
end1_index = iff.tindex_min(end_time1, timeB_us)


start_time2 = 60
end_time2 = 100
time_range2 = [start_time2,end_time2]

"""
start_time_index = iff.tindex_min(timeB_us[0,:],start_time)
end_time_index = iff.tindex_min(timeB_us[0,:],end_time)
"""
start_time_index = iff.tindex_min(start_time,timeB_us[0,:])
end_time_index = iff.tindex_min(end_time,timeB_us[0,:])

start_time_index2 = iff.tindex_min(start_time, timeB_us2[0,:])
end_time_index2 = iff.tindex_min(end_time, timeB_us2[0,:])

# Select shots to analyize
first_shot = 1
last_shot = 204
numshots = (last_shot-first_shot)+1
shot_range = [first_shot,last_shot]

first_shot2 = 1
last_shot2 = 145
numshots2 = (last_shot2-first_shot2)+1
shot_range2 = [first_shot2,last_shot2]

##np.zeros([THIS NUMBER,...]) needs to match poslist index
Bmag = np.zeros([6,numshots,len(data['pos5']['Bfilt']['z'][0,:])])
for pos in np.arange(len(poslist)):
        #x=data[poslist[pos]]['Bfilt']['r'][numshots,:]
        #y=data[poslist[pos]]['Bfilt']['z'][numshots,:]
        x=data[poslist[pos]]['Bfilt']['r'][:]
        y=data[poslist[pos]]['Bfilt']['z'][:]
        Bmag[pos,:] = np.sqrt(x**2+y**2)
        
Bmag2 = np.zeros([6,numshots2,len(data2['pos5']['Bfilt']['z'][0,:])])
for pos in np.arange(len(poslist2)):
        x2=data2[poslist2[pos]]['Bfilt']['r'][:]
        y2=data2[poslist2[pos]]['Bfilt']['z'][:]
        Bmag2[pos,:] = np.sqrt(x2**2+y2**2)        
        
fsize=int((data['pos5']['Bfilt']['z'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot = np.zeros([6,2,fsize])
avebspec_direct = np.zeros([6,2,fsize])
avebmagspec = np.zeros([6,fsize])

fsize2=int((data2['pos5']['Bfilt']['z'][0,start_time_index:end_time_index].shape[0])/2)+1
avebspec_frombdot2 = np.zeros([4,2,fsize2])
avebspec_direct2 = np.zeros([4,2,fsize2])
avebmagspec2 = np.zeros([4,fsize2])

dirlist = ['r','z']
dirlist2 = ['r','z']

for shot in np.arange(first_shot-1,last_shot):
    for pos in np.arange(len(poslist)):
        for direction in np.arange(len(dirlist)):
            """
            THE FIX FOR AMELIA IS CHANGING THE TIME_S[STARTTIMEINDEX:ENDTIMEINDEX] TO TIME_S[0,STARTTIME:ENDTIME]
            """
            """
            For fast velocities, lines below, look at spectra from index 68 and 72, plot fast and slow on one plot.
            """
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['Bdot'][dirlist[direction]][72,start_time_index:end_time_index],time_s[0,start_time_index:end_time_index],window='hanning')
            avebspec_frombdot[pos,direction,:]=avebspec_frombdot[pos,direction]+(pwr/(f*f))
        
            f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(data[poslist[pos]]['Bfilt'][dirlist[direction]][72,start_time_index:end_time_index],timeB_s[0,start_time_index:end_time_index],window='hanning')
            avebspec_direct[pos,direction,:]=avebspec_direct[pos,direction]+(pwr)
            
        
            
for shot in np.arange(first_shot2-1,last_shot2):
    for pos in np.arange(len(poslist2)):
        for direction in np.arange(len(dirlist2)):
            
            """
            For slow velocities, lines below, look at spectra from shots 135 and 135, plot fast and slow on one plot.
            """
            
            f2,f02,comp12,pwr2,mag12,phase12,cos_phase12,interval2=spec.spectrum_wwind(data2[poslist2[pos]]['Bdot'][dirlist2[direction]][135,start_time_index:end_time_index],time_s2[0,start_time_index:end_time_index],window='hanning')
            avebspec_frombdot2[pos,direction,:]=avebspec_frombdot2[pos,direction]+(pwr2/(f2*f2))
            f2,f02,comp12,pwr2,mag12,phase12,cos_phase12,interval2=spec.spectrum_wwind(data2[poslist2[pos]]['Bfilt'][dirlist2[direction]][135,start_time_index:end_time_index],timeB_s2[0,start_time_index:end_time_index],window='hanning')
            avebspec_direct2[pos,direction,:]=avebspec_direct2[pos,direction]+(pwr2)
            
            
#for shot in np.arange(first_shot-1,last_shot):
for shot in np.arange(numshots):
    for pos in np.arange(len(poslist)):
        f,f0,comp1,pwr,mag1,phase1,cos_phase1,interval=spec.spectrum_wwind(Bmag[pos,72,start_time_index:end_time_index],timeB_s[0,start_time_index:end_time_index],window='hanning')
        avebmagspec[pos,:]=avebmagspec[pos,:]+pwr

#for shot in np.arange(first_shot2-1,last_shot2):
for shot in np.arange(numshots2):
    for pos in np.arange(len(poslist2)):
        f2,f02,comp12,pwr2,mag12,phase12,cos_phase12,interval2=spec.spectrum_wwind(Bmag2[pos,135,start_time_index:end_time_index],timeB_s2[0,start_time_index:end_time_index],window='hanning')
        avebmagspec2[pos,:]=avebmagspec2[pos,:]+pwr2
        
                
allpos_avebspec_frombdot_2 = np.sum(avebspec_frombdot,axis=0)
allpos_avebspec_frombdot_2 = np.sum(allpos_avebspec_frombdot_2,axis=0)
allpos_avebmagspec_2 = np.sum(avebmagspec,axis=0)

allpos_avebspec_frombdot_22 = np.sum(avebspec_frombdot2,axis=0)
allpos_avebspec_frombdot_22 = np.sum(allpos_avebspec_frombdot_22,axis=0)
allpos_avebmagspec_22 = np.sum(avebmagspec2,axis=0)



###### Plot Details: #################
#plt.rc('axes',linewidth=2.0)
#plt.rc('xtick.major',width=2.0)
#plt.rc('ytick.major',width=2.0)
#plt.rc('xtick.minor',width=2.0)
#plt.rc('ytick.minor',width=2.0)
#plt.rc('lines',markersize=2.0,markeredgewidth=0.0,linewidth=1.0)
#fig=plt.figure(num=1,figsize=(4,3),dpi=300,facecolor='w',edgecolor='k')
#left  = 0.2  # the left side of the subplots of the figure
#right = 0.94    # the right side of the subplots of the figure
#bottom = 0.2  # the bottom of the subplots of the figure
#top = 0.90      # the top of the subplots of the figure
#wspace = 0.2   # the amount of width reserved for blank space between subplots
#hspace = 0.25   # the amount of height reserved for white space between subplots
#plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
#ax1=plt.subplot(1,1,1)   
######################################    
            
# Compute indices from freq range
start_f = 1e5
end_f = 7e5
start_f_index = iff.tindex_min(start_f,f)
end_f_index = iff.tindex_min(end_f,f)

start_f_index2 = iff.tindex_min(start_f,f2)
end_f_index2 = iff.tindex_min(end_f,f2)

dlog=np.log10(allpos_avebmagspec_2[start_f_index:end_f_index])
flog=np.log10(f[start_f_index:end_f_index])
A1=np.array([flog,np.ones(len(flog))])
w1=np.linalg.lstsq(A1.T,dlog)[0]
slope=np.round(w1[0],3)
func = f[1:]**(slope)
scale_data=allpos_avebmagspec_2[100]
scale_fit=func[100]
ratio=scale_data/scale_fit

dlog2=np.log10(allpos_avebmagspec_22[start_f_index2:end_f_index2])
flog2=np.log10(f2[start_f_index2:end_f_index2])
A12=np.array([flog2,np.ones(len(flog2))])
w12=np.linalg.lstsq(A12.T,dlog2)[0]
slope2=np.round(w12[0],3)
func2 = f2[1:]**(slope2)
scale_data2=allpos_avebmagspec_22[100]
scale_fit2=func2[100]
ratio2=scale_data2/scale_fit2

# I'm going to convert the frequency to a spatial scale by diving by speed v = 44km/s = 44000m/s
# This is the -5/3 Kolmogorov scaling we can compare our spectra to
#x_kolmo = f[1:]
#y_kolmo = np.power(f[1:], -5/3)

x_kolmo_freq = f[1:]
y_kolmo_freq = np.power(f[1:]/44000, -5/3)
x_kolmo_length = f[1:]/44000
y_kolmo_length = np.power(f[1:]/44000, -5/3)

Tmicroscale_temporal = np.zeros(6251)
Tmicroscale_spatial = np.zeros(6251)

Tmicroscale_temporal_freq = Tmicroscale_temporal+(3.51e-2**-1*44000) #Hertz
Tmicroscale_spatial_freq = Tmicroscale_spatial+(2.3e-2**-1*44000) #Hesrtz

Tmicroscale_temporal_length = Tmicroscale_temporal+(3.51e-2**-1) #inverse cm
Tmicroscale_spatial_length = Tmicroscale_spatial+(2.3e-2**-1) #inverse cm


#fig1, (ax1, ax2) = plt.subplots(2)
fig1, (ax1) = plt.subplots(1)
#plt.loglog(f2[1:],allpos_avebmagspec_22[1:],label='0.5ms')
#plt.loglog(f[start_f_index:end_f_index],allpos_avebmagspec_2[start_f_index:end_f_index])
#plt.loglog(f2[start_f_index2:end_f_index2],allpos_avebmagspec_22[start_f_index2:end_f_index2])

ax1.loglog(f[1:5000]/54000,allpos_avebmagspec_2[1:5000], color='blue', label='54 km/s', alpha=0.5)
ax1.loglog(f[1:5000]/211000, allpos_avebmagspec_22[1:5000], color='green', alpha=0.5, label='211 km/s') #label='211 km/s', alpha=0.5)
ax1.legend(loc='best', frameon=False)
#ax1.loglog(f[1:5000],allpos_avebspec_frombdot_2[1:5000], color='blue')
#ax1.axvline(30000, linestyle='--', color='black') #This is the 30kHz frequency filter (Low limit)
#ax1.axvline(13000000, linestyle=':', color='red') #This is the 13MHz Ion-Gyrofrequency (high limit)
#ax1.set_xlim(8000,60000000)
#ax1.set_ylim(0.0000000000001,10)
#ax1.loglog(x_kolmo_freq, y_kolmo_freq, color='black', linewidth=1)
#ax1.loglog(Tmicroscale_spatial_freq[3500:], allpos_avebmagspec_2[3500:], linestyle='-', color='black')
#ax1.loglog(Tmicroscale_temporal_freq[3500:], allpos_avebmagspec_2[3500:], linestyle='-', color='black')
ax1.set_xlabel(r'Frequency (Hz)')
ax1.set_ylabel(r'B-field Power (arb)')
#ax1.legend(loc='best', bbox_to_anchor=(0.25,0.3))
#ax1.text(2250000, 1e-12, r'$\lambda_{Taylor}^{Spatial}$')
#ax1.text(355000, 1e-12, r'$\lambda_{Taylor}^{Temporal}$')
#ax1.text(0.32, 0.95, 'Kolmogorov -5/3 Scale', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
#ax1.text(0.9, 0.9, '(a)', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
#plt.text(0.07,0.92,'(a)',horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontsize=12)
#plt.text(0.89, 0.00000000000001, r'$\frac{1}{m}$', fontsize=10)
#plt.text(8, 0.00000000000001, r'$\frac{1}{10cm}$', fontsize=10)
#plt.text(80, 0.00000000000001, r'$\frac{1}{cm}$', fontsize=10)
#plt.text(43, 0.000000000001, r'$\lambda <sub>t</sub>$')

####*****The two lines below add the orange line to the spectra plot!!! Removing f[start_f_index:end_f_index] plot removes the discoloration*****#####
#plt.loglog(f[start_f_index:end_f_index]/44000,allpos_avebmagspec_2[start_f_index:end_f_index])
#plt.loglog(f2[start_f_index2:end_f_index2]/44000,allpos_avebmagspec_22[start_f_index2:end_f_index2])

#ax2.loglog(f[1:5000]/44000,allpos_avebmagspec_2[1:5000], color='red')
#ax2.set_xlim(float(0.182), float(1363.636))
#plt.loglog(f2[1:]/44000,allpos_avebmagspec_22[1:],label='0.5ms')
#ax2.loglog(x_kolmo_length, y_kolmo_length, color='black', linewidth=1)

# The leftmost one is the temporal taylormicroscale
# The rightmost one is the spatial taylormicroscale

#ax2.loglog(Tmicroscale_spatial_length[3500:], allpos_avebmagspec_2[3500:], linestyle='-', linewidth=float(1.2), color='black', alpha=1)
#ax2.loglog(Tmicroscale_temporal_length[3500:], allpos_avebmagspec_2[3500:], linestyle='-', linewidth=float(1.2), color='black', alpha=1)
#plt.xticks(fontsize=9)
#plt.yticks(fontsize=9)
#ax2.set_xlabel(r'Inverse Length ($m^{-1})$')
#ax2.set_ylabel('B-field Power (arb)')
#plt.grid(True, alpha=0.5)
#ax2.text(0.7, 4, 'Kolmogorov -5/3 Scale')
#ax2.text(0.32, 0.95, 'Kolmogorov -5/3 Scale', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
#ax2.text(0.9, 0.9, '(b)', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
#ax2.text(4, 2e-15, r'$\frac{1}{25cm}$')
#ax2.text(0.895, 2e-15, r'$\frac{1}{m}$')
#ax2.text(8, 2e-15, r'$\frac{1}{10cm}$')
#ax2.text(80, 2e-15, r'$\frac{1}{cm}$')
#ax2.text(50, 1e-12, r'$\lambda_{Taylor}^{Spatial}$')
#ax2.text(7.5, 1e-12, r'$\lambda_{Taylor}^{Temporal}$')
#ax2.legend(loc='best', bbox_to_anchor=(0.25,0.3))
#plt.legend(loc='best')


print ('1.5ms slope is',slope)
print ('0.5ms slope is',slope2)

#plt.savefig('C:\\Users\\Josh0\\Documents\\1. Josh Documents\\Graduate School - Bryn Mawr College\\Plasma Lab (BMX) Research\\Analysis\\Plots\\9222021_Spectra_1p5delay_30khz_13MHz.pdf', dpi=600)
            
            
            
            
            
            
            
            
            