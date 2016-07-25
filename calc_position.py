from __future__ import division
import numpy as np
import os
from astropy.time import Time
import ephem

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



