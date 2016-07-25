from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import gridspec
from functions_idl_image import *
from astropy.convolution import convolve, Box1DKernel


"""
The following code allows a user to set a time and thus calculate the positions of the satellites
and Venus. The positions can then be used to generate a los-int image from IDL which outputs data 
to python and then an image of the disk at that time is obtained, along with a polar plot showing 
which part of the ring is being surveyed by the satellites for the given time
"""

time = '2009-06-07'	#date on which the circumsolar ring is observed, in yyyy/mm/dd format, dd can be decimal
sat = 'B' 	#takes the argument A or B

calc_ephem(time)

#a plot to show where Venus and the satellites are at the given time
#and the tangent point one is looking at through the satellite
tangent_a = stereo_a.hlon-np.arctan(0.72/np.sqrt((stereo_a.sun_distance**2)-(0.72**2)))
tangent_b = stereo_b.hlon+np.arctan(0.72/np.sqrt((stereo_b.sun_distance**2)-(0.72**2)))

if (sat == 'A'):
	filein = 'volimage.dat'
	fileout = 'output.sav'
	a = stereo_a.sun_distance
	e = STEREO_A.ecc
	i = STEREO_A.inc
	ilng = -60	
	flng = -35
	asc = princ_range(STEREO_A.asc-(ven.hlon*(180/np.pi)))		
	wbar = princ_range(STEREO_A.wbar-(ven.hlon*(180/np.pi))) 	
	inter_long = (stereo_a.hlon-ven.hlon)*(180/np.pi)
	meanlong = princ_range(inter_long)			#rel angle: venus-sun-sat angle
	long1 = inter_long + 180 + ilng
	long2 = inter_long + 180 + flng

elif (sat == 'B'):
	filein = 'volimage.dat'
	fileout = 'output.sav'
	a = stereo_b.sun_distance
	e = STEREO_B.ecc
	i = STEREO_B.inc
	ilng = 35
	flng = 60
	asc = princ_range(STEREO_B.asc-(ven.hlon*(180/np.pi)))
	wbar = princ_range(STEREO_B.wbar-(ven.hlon*(180/np.pi))) 
	inter_long = (stereo_b.hlon-ven.hlon)*(180/np.pi)
	meanlong = princ_range(inter_long)			#rel angle: venus-sun-sat angle
	long1 = inter_long + 180 + ilng
	long2 = inter_long + 180 + flng

output = losint_image(filein,fileout,a,e,i,asc,wbar,meanlong,lon1=long1,lon2=long2)

image = output[0]
lg = output[1]
lat = output[2]

filter_width = 13      #assign the filter width in pixels/cells
model_filtered = np.zeros((len(image),len(image[0])))
for i in range(0,len(image)):
	model_filtered[i] = convolve(image[i,:],Box1DKernel(filter_width))

fig = plt.figure(figsize=(15,8))
gs = gridspec.GridSpec(5,3)
ax1 = fig.add_subplot(gs[0,0])
ax1.plot(np.linspace(ilng,flng,50),np.average(image,0))
plt.locator_params(axis='y',nbins=5)
plt.setp(ax1.get_xticklabels(),visible=False)
ax1.set_title('Raw Image')
ax2 = fig.add_subplot(gs[1:5,0],sharex=ax1)
im = ax2.imshow(image,origin='lower',extent=(ilng,flng,lat[0],lat[-1]),aspect='auto')
ax2.set_xlabel('Helioecliptic Longitude')
ax2.set_ylabel('Helioecliptic Latitude')
plt.colorbar(im,orientation='horizontal')

ax3 = fig.add_subplot(gs[0,1])
ax3.plot(np.linspace(ilng,flng,50),np.average(model_filtered,0))
plt.locator_params(axis='y',nbins=5)
plt.setp(ax3.get_xticklabels(),visible=False)
ax3.set_title('Boxcar Filter')
ax4 = fig.add_subplot(gs[1:5,1],sharex=ax3)
im2 = ax4.imshow(model_filtered,origin='lower',extent=(ilng,flng,lat[0],lat[-1]),aspect='auto')
ax4.set_xlabel('Helioecliptic Longitude')
ax4.set_ylabel('Helioecliptic Latitude')
plt.colorbar(im2,orientation='horizontal')

final_image = np.subtract(image,model_filtered)
ax5 = fig.add_subplot(gs[0,2])
ax5.plot(np.linspace(ilng+3,flng-3,38),np.average(final_image[:,6:44],0))
plt.locator_params(axis='y',nbins=5)
plt.setp(ax5.get_xticklabels(),visible=False)
ax5.set_title('Final Image')
ax6 = fig.add_subplot(gs[1:5,2],sharex=ax5)
im3 = ax6.imshow(final_image[:,6:44],origin='lower',extent=(ilng+3,flng-3,lat[0],lat[-1]),aspect='auto')
ax6.set_xlabel('Helioecliptic Longitude')
ax6.set_ylabel('Helioecliptic Latitude')
plt.colorbar(im3,orientation='horizontal')

os.chdir('/home/yc367/Desktop/Test/Output')
plt.savefig(time+'_'+sat+'.png')

plt.figure(figsize=(7,4))
venus, = plt.polar(ven.hlon,ven.sun_distance,'o',color='y')
a, = plt.polar(stereo_a.hlon,stereo_a.sun_distance,'o',color='b')
tan_a, = plt.polar(tangent_a,0.72,'s',color='b')
b, = plt.polar(stereo_b.hlon,stereo_b.sun_distance,'o',color='r')
tan_b, = plt.polar(tangent_b,0.72,'s',color='r')
plt.polar(np.linspace(0,2*np.pi),[0.72]*50,linestyle='dashed',color='y')
plt.legend([venus,a,tan_a,b,tan_b],['Venus','STEREO A','a','STEREO B','b'],loc=(1.1,0.35),numpoints=1)
plt.title(time)
plt.ylim([0,1.2])
plt.tight_layout()
plt.savefig(time+'_'+'Positions.png')
