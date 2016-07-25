from __future__ import division
import numpy as np
import os
from astropy.time import Time
import ephem
from subprocess import Popen,PIPE
from scipy.io import readsav
from os import remove
from os.path import isfile

def princ_range(angle):
	if angle < 0:
		return angle+360
	if angle > 360:
		return angle-360
	else:
		return angle

class STEREO_B: 
	a = 1.0428910249590357
	ecc = 0.041433775958706259
	inc = 0.29399937403029891
	asc = 336.41463527732202
	peri = 146.96727037160804
	wbar = princ_range(asc+peri)

	b_orb = 'STEREO B,e,'+str(inc)+','+str(asc)+','+str(peri)+','+str(a)+',0.925759736,'+str(ecc)+',186.385347,08/13/2007,2000,g 10.0,1.0'

stereo_b = ephem.readdb(STEREO_B.b_orb)

class STEREO_A:
	a = 0.96177266435859354
	ecc = 0.0057597588012127257
	inc = 0.12649088953666768
	asc = 214.00561578507018
	peri = 93.60768134897603
	wbar = princ_range(asc+peri)

	a_orb = 'STEREO A,e,'+str(inc)+','+str(asc)+','+str(peri)+','+str(a)+',1.044486455,'+str(ecc)+',26.5232145634,08/13/2007,2000,g 10.0,1.0'

stereo_a = ephem.readdb(STEREO_A.a_orb)

ven = ephem.Venus()
earth = ephem.Sun()

def calc_ephem(time):
	ven.compute(time)
	earth.compute(time)
	stereo_a.compute(time)
	stereo_b.compute(time)


'''Call IDL to generate a save file from the given ring model and
    viewing geometry'''
def losint_image(filein,fileout,a,e,i,asc,wbar,meanlong,
                 nlon=50,nlat=17,lon1=35.,lon2=60.,lat1=-8.5,lat2=8.5,nstep=200):

    if isfile(fileout):
        print("ERROR: Output file exists already, exiting.")
        return

    # run IDL to generate a save file
    args = 'ring_image,"'+filein+'","'+fileout+'",'+str(a)+','+str(e)+','+str(i)+','+str(asc)+','+str(wbar)+','+str(meanlong)+',lon1='+str(lon1)+',lon2='+str(lon2)+',lat1='+str(lat1)+',lat2='+str(lat2)+',nlon='+str(nlon)+',nlat='+str(nlat)
    p1 = Popen(['echo',args],stdout=PIPE)
    p2 = Popen('idl',stdin=p1.stdout,stdout=PIPE)
    p1.stdout.close()
    aout = p2.communicate()
    print(aout)

    # grab desired output from save file and then delete it
    out = readsav(fileout)
    remove(fileout)
    return (out.im,out.lons,out.lats)
