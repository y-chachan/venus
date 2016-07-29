
#!/bin/env python
import pickle
import numpy as np
import matplotlib.pyplot as plt
import os,glob
from matplotlib import gridspec
import operator
from astropy.convolution import convolve, Box1DKernel
from matplotlib import ticker
from functions import *

#a piece of code to make cumulative plots of means scans, final images, etc
path = 'C:/Users/YaChan/Desktop/Output_A'   #place where input is, the output pickle files
sat = 'A'
os.chdir(path)
fnames = glob.glob('*.p')

arr = np.zeros((len(fnames),17,50))
for i in range(0,len(fnames)):
    arr[i] = pickle.load(open(fnames[i],'rb'))

fig = plt.figure(figsize=(5,len(fnames)))
for i in range(0,len(fnames)):    
    mean_scan = np.average(arr[i],0)
    if sat == 'A':
        lg = np.linspace(-60,-35,50)
        plt.xlim([-57,-38])
    if sat == 'B':
        lg = np.linspace(35,60,50)
        plt.xlim([38,57])
    plt.plot(lg[6:44],mean_scan[6:44]+i/40)
    plt.hlines(i/40,lg[6],lg[44],linestyles='dashed',color='0.75')
    plt.annotate(fnames[i].split('.')[0],xy=(0.8,(i+1)/(len(fnames)+1)),xycoords='axes fraction',size='small')
    
plt.ylim([-0.025,len(fnames)/40])
plt.xlabel('Helioecliptic Longitude')
os.chdir('C:/Users/YaChan/Desktop')
plt.savefig('Mean_Scans_'+sat+'.png')

fig2 = plt.figure(figsize=(6,len(fnames)*4))
gs2 = gridspec.GridSpec(len(fnames),1)
for i in range(0,len(fnames)):      
    ax = fig2.add_subplot(gs2[i,0])
    if sat == 'A':
        extent = (-57,-38,-8.5,8.5)
    if sat == 'B':
        extent = (38,57,-8.5,8.5)
    im = ax.imshow(arr[i][:,6:44],origin='lower',extent=extent,aspect='auto')
    plt.setp(ax.get_xticklabels(),visible=False)
    plt.title(fnames[i].split('.')[0])
    plt.colorbar(im)
#plt.tight_layout()    
plt.setp(ax.get_xticklabels(),visible=True)
plt.xlabel('Helioecliptic Longitude')
os.chdir('C:/Users/YaChan/Desktop')
plt.savefig('Images_'+sat+'.png')

#piece of code to find midrise points and mid-drop points, plot them against distance
#can also be used to see how the peak height and width evolve with time
path = 'C:/Users/YaChan/Desktop/Out/B'   #place where input is, the output pickle files
os.chdir(path)
sat = 'B'
fnames = glob.glob('*.p')
file_count = len(fnames)
date = np.zeros((file_count)).astype('str')
index_midrise = np.zeros((file_count)).astype('int')
index_midfall = np.zeros((file_count)).astype('int')
sat_sun_distance = np.zeros((file_count))


for i in range(0,file_count):
    input  = pickle.load(open(fnames[i], 'rb'))
    date[i] = fnames[i].split('.')[0]
    calc_ephem(date[i])
    if (sat=='A'):
        sat_sun_distance[i] = stereo_a.sun_distance
    elif (sat=='B'):
        sat_sun_distance[i] = stereo_b.sun_distance
    
    final_image = input['final_image']
    #necessary to truncate because we are trying to find the midrise and mid-drop points,
    #need to identify maxima and minima correctly for that
    sub_image = final_image[:,6:44] 
    mean_scan = np.average(sub_image,0)

    index_max, value_max = max(enumerate(mean_scan),key=operator.itemgetter(1))     #can store value_max to store peak heights
    #print(index_max,lg[index_max],value_max)
    index_min1, value_min1 = min(enumerate(mean_scan[:index_max]),key=operator.itemgetter(1))
    #print(index_min1,lg[index_min1],value_min1)
    index_min2, value_min2 = min(enumerate(mean_scan[index_max:]),key=operator.itemgetter(1))
    #print(index_min2+index_max,lg[index_min2+index_max],value_min2)
    if (sat=='B'):
        midrise = (value_min1+value_max)/2
        midfall = (value_max+value_min2)/2

        index_midrise[i] = min(range(index_min1,index_max), key=lambda i: abs(mean_scan[i]-midrise))
        index_midfall[i] = min(range(index_max,index_max+index_min2), key=lambda i: abs(mean_scan[i]-midfall))
    if (sat=='A'):
        midrise = (value_min2+value_max)/2
        midfall = (value_max+value_min1)/2

        index_midrise[i] = min(range(index_max,index_max+index_min2), key=lambda i: abs(mean_scan[i]-midrise))
        index_midfall[i] = min(range(index_min1,index_max), key=lambda i: abs(mean_scan[i]-midfall))

    #code to plot the midrise and midfall points on the plot 
    """plt.figure(figsize=(5,2))
    plt.plot(lg,mean_scan)
    plt.axhline(linestyle='dashed',color='0.75')
    plt.axvline(lg[index_midrise],0,1,linestyle='dashed',color='0.75')
    plt.axvline(lg[index_midfall],0,1,linestyle='dashed',color='0.75')
    plt.title('Mean Scan')
    plt.ylabel('DN s$^{-1}$pix$^{-1}$')
    plt.xlabel('Helioecliptic Longitude')
    plt.xlim([38,57])"""

#peak_widths = np.subtract(lg[index_midfall],lg[index_midrise])
#a plot of the longtiude of midrise against the distance of satellite at a given time
if (sat=='A'):
    lg = np.linspace(-57,-38,38)  #these will change depending on the filter width of boxcar and the extent of the image
    h_dist = np.linspace(STEREO_A.a*(1-STEREO_A.ecc),STEREO_A.a*(1+STEREO_A.ecc),1000)
elif (sat=='B'):
    lg = np.linspace(38,57,38)  #these will change depending on the filter width of boxcar and the extent of the image
    h_dist = np.linspace(STEREO_B.a*(1-STEREO_B.ecc),STEREO_B.a*(1+STEREO_B.ecc),1000)

plt.figure()
inv = 1./h_dist
graph1 = np.rad2deg(np.arcsin(0.70*inv))
graph2 = np.rad2deg(np.arcsin(0.72*inv))
graph3 = np.rad2deg(np.arcsin(0.71*inv))
graph4 = np.rad2deg(np.arcsin(0.73*inv))
plt.plot(h_dist,graph1)
plt.plot(h_dist,graph2,linestyle='dashed')
plt.plot(h_dist,graph3,linestyle='dashed')
plt.plot(h_dist,graph4)
#plt.annotate('0.7 AU',xy=(0.85,0.05),xycoords='axes fraction',size='small')
#plt.annotate('0.73 AU',xy=(0.9,0.35),xycoords='axes fraction',size='small')
plt.scatter(sat_sun_distance,abs(lg[index_midrise])) #,yerr=0.5,fmt='o')
plt.xlabel('Heliocentric Distance (AU)')
plt.ylabel('Helioecliptic Longitude (degrees)')
plt.xlim([h_dist[0],h_dist[-1]])
plt.title('STEREO '+sat)
plt.savefig('Long vs Dist_'+sat)