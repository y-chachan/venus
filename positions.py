from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as mpatches
from calc_position import *

#path2 = '/home/yc367/Desktop/Positions_ephem'  
#os.chdir(path2)

time = np.arange('2007-08-01', '2016-06-30', dtype='datetime64[W]')

milky_way_a = np.zeros((10,2)).astype(object)
milky_way_a[:,0] = np.arange(np.datetime64('2007-09-04'),np.datetime64('2016-06-30'),np.timedelta64(346,'D'))
milky_way_a[:,1] = np.arange(np.datetime64('2007-11-04'),np.datetime64('2016-06-30'),np.timedelta64(346,'D'))

milky_way_b = np.zeros((10,2)).astype(object)
milky_way_b[:,0] = np.arange(np.datetime64('2008-01-31'),np.datetime64('2018-06-30'),np.timedelta64(388,'D'))
milky_way_b[:,1] = np.arange(np.datetime64('2008-04-01'),np.datetime64('2018-06-30'),np.timedelta64(388,'D'))

hae_ven = np.zeros((len(time)))    
hae_earth = np.zeros((len(time)))
hEcl_Lon_stereo_a = np.zeros((len(time)))
hEcl_Lon_stereo_b = np.zeros((len(time))) 
ven_a = np.zeros((len(time)))
earth_a = np.zeros((len(time)))
stereoa_a = np.zeros((len(time)))
stereob_a = np.zeros((len(time)))

for i in range(0,len(time)):
    calc_ephem(str(time[i]))
    hEcl_Lon_stereo_a[i] = stereo_a.hlon
    hEcl_Lon_stereo_b[i] = stereo_b.hlon
    hae_ven[i] = ven.hlon
    hae_earth[i] = earth.hlon
    ven_a[i] = ven.sun_distance
    earth_a[i] = earth.earth_distance
    stereoa_a[i] = stereo_a.sun_distance
    stereob_a[i] = stereo_b.sun_distance


ven_elongation_a = np.subtract(hEcl_Lon_stereo_a,hae_ven)
earth_elongation_a = np.subtract(hae_earth,hEcl_Lon_stereo_a)
ven_elongation_b = np.subtract(hEcl_Lon_stereo_b,hae_ven)
earth_elongation_b = np.subtract(hae_earth,hEcl_Lon_stereo_b)

for i in range(0,len(time)):
    ven_elongation_a[i] = princ_range(ven_elongation_a[i],'r')
    earth_elongation_a[i] = princ_range(earth_elongation_a[i],'r')
    ven_elongation_b[i] = princ_range(ven_elongation_b[i],'r')
    earth_elongation_b[i] = princ_range(earth_elongation_b[i],'r')


ven_dist_a = np.sqrt(np.power(stereoa_a,2) + np.power(ven_a,2) - 2*np.multiply(np.multiply(stereoa_a,ven_a),np.cos(ven_elongation_a)))
ven_dist_b = np.sqrt(np.power(stereob_a,2) + np.power(ven_a,2) - 2*np.multiply(np.multiply(stereob_a,ven_a),np.cos(ven_elongation_b)))

ven_subtend_a = np.rad2deg(np.arcsin(np.divide(np.multiply(ven_a,np.sin(ven_elongation_a)),ven_dist_a)))
ven_subtend_b = np.rad2deg(np.arcsin(np.divide(np.multiply(ven_a,np.sin(ven_elongation_b)),ven_dist_b)))

earth_subtend_a = 0.5 * np.subtract(180,np.absolute(np.rad2deg(earth_elongation_a)))
earth_subtend_b = 0.5 * np.subtract(180,np.absolute(np.rad2deg(earth_elongation_b)))


tangent_a = np.subtract(hEcl_Lon_stereo_a,np.arctan(0.72/np.sqrt((np.power(stereoa_a,2)-(0.72**2)))))
tangent_b = np.add(hEcl_Lon_stereo_b,np.arctan(0.72/np.sqrt((np.power(stereob_a,2)-(0.72**2)))))
for i in range(0,len(time)):
    tangent_a[i] = princ_range(tangent_a[i],'r')
    tangent_b[i] = princ_range(tangent_b[i],'r')

tan_rel_a = np.subtract(tangent_a,hae_ven)
tan_rel_b = np.subtract(tangent_b,hae_ven)
for i in range(0,len(time)):
    tan_rel_a[i] = princ_range(tan_rel_a[i],'r')
    tan_rel_b[i] = princ_range(tan_rel_b[i],'r')

fig = plt.figure(figsize=(15,15))
gs = gridspec.GridSpec(4,1)
ax1 = fig.add_subplot(gs[0,0])
earth, = ax1.plot(time,earth_subtend_b)
ven, = ax1.plot(time,ven_subtend_b)
ax1.hlines([40,50],time[0],time[-1],linestyles='dashed',color = '0.75')
for i in range(0,10):
    ax1.axvspan(milky_way_b[i,0],milky_way_b[i,1], alpha=0.25, color='k')
ax3 = fig.add_subplot(gs[1,0],sharex=ax1)
ax3.plot(time,np.rad2deg(tan_rel_b))
#ax3.plot(time,np.rad2deg(ven_elongation_b))
ax2 = fig.add_subplot(gs[2,0],sharex=ax1)
ax2.plot(time,earth_subtend_a)
ax2.plot(time,ven_subtend_a)
ax2.hlines([-50,-40],time[0],time[-1],linestyles='dashed',color = '0.75')
for i in range(0,10):
    ax2.axvspan(milky_way_a[i,0],milky_way_a[i,1], alpha=0.25, color='k')
ax4 = fig.add_subplot(gs[3,0],sharex=ax1)
ax4.plot(time,np.rad2deg(tan_rel_a))
#ax4.plot(time,np.rad2deg(ven_elongation_a))
ax1.set_title('STEREO B')
ax2.set_title('STEREO A')
ax4.set_xlabel('Time')
ax3.set_xlabel('Time')
ax1.set_ylabel('Helioecliptic Longitude')
ax2.set_ylabel('Helioecliptic Longitude')
ax3.set_ylabel('Elongation of Tangent Point')
ax4.set_ylabel('Elongation of Tangent Point')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax1.set_xlim([time[0],time[-1]])
ax2.set_ylim([-90,0])
ax1.set_ylim([0,90])
ax3.set_ylim(0,360)
ax4.set_ylim(0,360)
ax3.grid()
ax4.grid()
minorLocator = AutoMinorLocator(n=12)
ax1.get_xaxis().set_minor_locator(minorLocator)
plt.tight_layout()
grey_patch = mpatches.Patch(color='0.75', label='Milky Way')
plt.legend([ven, earth, grey_patch], ['Venus', 'Earth', 'Milky Way'],bbox_to_anchor=(0.65, 4.75, 1., .102), loc=3,
           ncol=3, borderaxespad=0.)
plt.savefig('HPLN- Venus & Earth.png')


fig2 = plt.figure(figsize=(25,4))
gs = gridspec.GridSpec(1,7)
ax1 = fig2.add_subplot(gs[0,0:3])
#ax1.plot(time,np.rad2deg(ven_elongation_b))
ax1.plot(time,np.rad2deg(tan_rel_b))
ax1.axvspan('2008-04-01','2008-09-10', alpha=0.25, color='r')	#STEREO B
ax1.axvspan('2009-06-07','2009-08-07', alpha=0.25, color='r')	#STEREO B
ax1.axvspan('2014-02-01','2014-04-30', alpha=0.25, color='r')	#STEREO B
ax1.axhspan(130,230, alpha=0.25, color='r')
ax1.axhspan(260,300, alpha=0.25, color='r')
ax2 = fig2.add_subplot(gs[0,3:6],sharey=ax1)
#ax2.plot(time,np.rad2deg(ven_elongation_a))
ax2.plot(time,np.rad2deg(tan_rel_a))
ax1.set_title('STEREO B')
ax2.set_title('STEREO A')
ax1.set_ylabel('Venus Elongation')
ax1.set_xlabel('Time')
ax2.set_xlabel('Time')
ax1.set_ylim(0,360)
ax2.set_ylim(0,360)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.axvspan('2009-06-07','2009-08-07', alpha=0.25, color='b')	#STEREO A
#ax2.axvspan(2009+(158/365),2009+(199/365), alpha=0.25, color='b')	#STEREO A
ax2.axvspan('2013-02-01','2013-04-30', alpha=0.25, color='b')	#STEREO A
ax2.axhspan(280,320, alpha=0.25, color='b')
ax2.axhspan(240,280, alpha=0.25, color='b')
ax3 = fig2.add_subplot(gs[0,6],sharey=ax1)
ax3.axhspan(130,230, alpha=0.25, color='r')
ax3.axhspan(260,300, alpha=0.25, color='r')
ax3.axhspan(280,320, alpha=0.25, color='b')
ax3.axhspan(240,280, alpha=0.25, color='b')
ax3.set_ylim(0,360)
plt.savefig('Analysed')
"""
plt.xticks(x,tx,rotation='vertical')
plt.annotate('Venus Circumsolar Ring',xy=(0.65,0.71),xycoords='axes fraction')
plt.annotate('Venus Circumsolar Ring',xy=(0.65,0.26),xycoords='axes fraction')"""
