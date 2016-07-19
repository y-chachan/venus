from __future__ import division     #only necessary for older python versions
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.time import Time
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as mpatches
import ephem


path2 = '/home/yc367/Desktop/Positions_ephem'  #set the path to where the txt files are
os.chdir(path2)

stereo_a = np.genfromtxt('stereoa.txt',skip_header=50,skip_footer=32)
stereo_b = np.genfromtxt('stereob.txt',skip_header=50,skip_footer=32)

time = stereo_a[:,0]    #julian days
t = Time(time,format='jd')
time_dy = t.decimalyear

dublin_start_date = ephem.julian_date('1899/12/31 12:00:00')
time_ephem = np.subtract(time,dublin_start_date)    #ephem counts dates from the dublin start date

hEcl_Lon_stereo_a = stereo_a[:,1]
hEcl_Lon_stereo_b = stereo_b[:,1]       

milky_way_a = np.zeros((10,2))
milky_way_a[0] = (2007+(247/365),2007+(308/365))
for i in range(1,10):
    milky_way_a[i] = np.add(milky_way_a[i-1],346/365)

milky_way_b = np.zeros((10,2))
milky_way_b[0] = (2008+(31/365),2008+(92/365))
for i in range(1,10):
    milky_way_b[i] = np.add(milky_way_b[i-1],388/365)

ven = ephem.Venus()
earth = ephem.Sun()

hae_ven = np.zeros((len(time)))    
hae_earth = np.zeros((len(time)))

for i in range(0,len(time)):
	ven.compute(time_ephem[i])
	earth.compute(time_ephem[i])
	hae_ven[i] = ven.hlon*(180/np.pi)
	hae_earth[i] = earth.hlon*(180/np.pi)

ven_elongation_a = np.subtract(hae_ven,hEcl_Lon_stereo_a)
earth_elongation_a = np.subtract(hae_earth,hEcl_Lon_stereo_a)

ven_elongation_b = np.subtract(hae_ven,hEcl_Lon_stereo_b)
earth_elongation_b = np.subtract(hae_earth,hEcl_Lon_stereo_b)

for i in range(0,len(time)):
    if ven_elongation_a[i] < 0:  ven_elongation_a[i] += 360
    if earth_elongation_a[i] < 0:  earth_elongation_a[i] += 360
    if ven_elongation_b[i] < 0:  ven_elongation_b[i] += 360
    if earth_elongation_b[i] < 0:  earth_elongation_b[i] += 360
        
ven_dist_a = np.sqrt(1 + (0.72**2) - 2*0.72*np.cos(ven_elongation_a*np.pi/180))
ven_dist_b = np.sqrt(1 + (0.72**2) - 2*0.72*np.cos(ven_elongation_b*np.pi/180))
ven_subtend_a = np.rad2deg(np.arcsin(0.72*np.sin(ven_elongation_a*np.pi/180)/ven_dist_a))
ven_subtend_b = np.rad2deg(np.arcsin(0.72*np.sin(ven_elongation_b*np.pi/180)/ven_dist_b))

earth_subtend_a = 0.5 * np.subtract(180,np.absolute(earth_elongation_a))
earth_subtend_b = 0.5 * np.subtract(180,np.absolute(earth_elongation_b))

fig = plt.figure(figsize=(15,15))
gs = gridspec.GridSpec(4,1)
ax1 = fig.add_subplot(gs[0,0])
earth, = ax1.plot(time_dy,earth_subtend_b)
ven, = ax1.plot(time_dy,ven_subtend_b)
ax1.hlines([40,50],time_dy[0],time_dy[-1],linestyles='dashed',color = '0.75')
for i in range(0,10):
    ax1.axvspan(milky_way_b[i,0],milky_way_b[i,1], alpha=0.25, color='k')
ax3 = fig.add_subplot(gs[1,0],sharex=ax1)
ax3.plot(time_dy,ven_elongation_b)
ax2 = fig.add_subplot(gs[2,0],sharex=ax1)
ax2.plot(time_dy,earth_subtend_a)
ax2.plot(time_dy,ven_subtend_a)
ax2.hlines([-50,-40],time_dy[0],time_dy[-1],linestyles='dashed',color = '0.75')
for i in range(0,10):
    ax2.axvspan(milky_way_a[i,0],milky_way_a[i,1], alpha=0.25, color='k')
ax4 = fig.add_subplot(gs[3,0],sharex=ax1)
ax4.plot(time_dy,ven_elongation_a)
ax1.set_title('STEREO B')
ax2.set_title('STEREO A')
ax4.set_xlabel('Time')
ax3.set_xlabel('Time')
ax1.set_ylabel('Helioecliptic Longitude')
ax2.set_ylabel('Helioecliptic Longitude')
ax3.set_ylabel('Venus Elongation')
ax4.set_ylabel('Venus Elongation')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax1.set_xlim([time_dy[0],time_dy[-1]])
ax2.set_ylim([-90,0])
ax1.set_ylim([0,90])
ax3.set_ylim(0,360)
ax4.set_ylim(0,360)
ax3.grid()
ax4.grid()
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
minorLocator = AutoMinorLocator(n=12)
ax1.get_xaxis().set_minor_locator(minorLocator)
plt.tight_layout()
grey_patch = mpatches.Patch(color='0.75', label='Milky Way')
plt.legend([ven, earth, grey_patch], ['Venus', 'Earth', 'Milky Way'],bbox_to_anchor=(0.65, 4.75, 1., .102), loc=3,
           ncol=3, borderaxespad=0.)
plt.savefig('HPLN- Venus & Earth.png')