from __future__ import division
import numpy as np
import os,glob
from functions import *
import functions
import pickle

"""
The following code allows a user to set a time and thus calculate the positions of the satellites
and Venus. The positions can then be used to generate a los-int image from IDL which outputs data 
to python and then an image of the disk at that time is obtained, along with a polar plot showing 
which part of the ring is being surveyed by the satellites for the given time
"""

sat = 'A' 	#takes the argument A or B
beta = '0.05' 	#the ratio of radiation forces to gravitational force, a measure of particle size
time = '2009-06-12'	#date on which the circumsolar ring is observed, in yyyy/mm/dd format, dd can be decimal

calc_ephem(time)

#a plot to show where Venus and the satellites are at the given time
#and the tangent point one is looking at through the satellite
tangent_a = stereo_a.hlon-np.arctan(0.72/np.sqrt((stereo_a.sun_distance**2)-(0.72**2)))
tangent_b = stereo_b.hlon+np.arctan(0.72/np.sqrt((stereo_b.sun_distance**2)-(0.72**2)))

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

if (sat == 'A'):
	filein = 'volimage.dat'
	fileout = 'output.sav'
	a = STEREO_A.a
	e = STEREO_A.ecc
	i = STEREO_A.inc
	ilng = -180	
	flng = 180
	asc = princ_range(STEREO_A.asc-(ven.hlon*(180/np.pi)))		
	wbar = princ_range(STEREO_A.wbar-(ven.hlon*(180/np.pi))) 	
	true_long = (stereo_a.hlon-ven.hlon)*(180/np.pi)
	meanlong = calc_meanlong(true_long,wbar,e)	
	long1 = true_long + 180 + ilng
	long2 = true_long + 180 + flng

elif (sat == 'B'):
	filein = 'volimage.dat'
	fileout = 'output.sav'
	a = stereo_b.sun_distance
	e = STEREO_B.ecc
	i = STEREO_B.inc
	ilng = -180
	flng = 180
	asc = princ_range(STEREO_B.asc-(ven.hlon*(180/np.pi)))
	wbar = princ_range(STEREO_B.wbar-(ven.hlon*(180/np.pi))) 
	true_long = (stereo_b.hlon-ven.hlon)*(180/np.pi)
	meanlong = calc_meanlong(true_long,wbar,e)				
	long1 = true_long + 180 + ilng
	long2 = true_long + 180 + flng

output = losint_image(filein,fileout,a,e,i,asc,wbar,meanlong,lon1=long1,lon2=long2)

image = output[0]
lng = output[1]
lat = output[2]
lg = np.linspace(ilng,flng,50)
vgridsize = len(image)
hgridsize = len(image[0])

data_files = pickle.load(open(time+'.p','rb'))
IV = data_files['pl']

pl_added = np.zeros((vgridsize,hgridsize))

for i in range(0,len(image)):
    pl = (1/float(beta))*IV[i,0]*np.power(np.absolute(lg),IV[i,1]) 	#need to change the premultiplication factor for different beta (either scale it up or create each case for beta
    pl_added[i] = np.add(image[i],pl)

output = process_image(pl_added,vgridsize,hgridsize,lg)

inter_pl_added = np.subtract(pl_added,output[0])
inter_pl_image = np.subtract(output[1],output[2])
final_image2 = np.subtract(inter_pl_added,inter_pl_image)

output_file = time+'_'+sat+'_'+beta
plot_output(final_image,vgridsize,lg,ilng,flng,lat[0],lat[-1],output_file)
plot_diagnostics(pl_added,output[0],output[1],output[2],inter_pl_added,inter_pl_image,final_image,vgridsize,lg,ilng,flng,output[3],output_file)
