from __future__ import division 	#to be used with older python versions
from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import optimize
import functions

os.chdir('/home/yc367/Desktop/Disk_STEREO_images/0.1')

ilng = -70
flng = 70
ilat = -40
flat = 40

image = readsav('res_disk_sv.sav')
arr = image.im
plt.figure()
plt.imshow(arr,origin='lower',extent=(ilng,flng,ilat,flat),aspect='auto')
plt.vlines([-57,-38,38,57],-8.5,8.5,linestyles='solid',colors='w')
plt.hlines([-8.5,8.5],-57,-38,linestyles='solid',colors='w')
plt.hlines([-8.5,8.5],38,57,linestyles='solid',colors='w')
plt.annotate('STEREO B',xy = (0.75,0.65),xycoords='axes fraction',size='small',color='w')
plt.annotate('STEREO A',xy = (0.09,0.65),xycoords='axes fraction',size='small',color='w')
plt.xlabel('Helioecliptic Longitude')
plt.ylabel('Helioecliptic Latitude')
plt.colorbar()
plt.savefig('res_disk_sv.png')

midlatpx = len(arr[:,0])/2
midlongpx = len(arr[0])/2

lat_px_deg = len(arr[:,0])/(flat-ilat)
long_px_deg = len(arr[0])/(flng-ilng)

arr_stereob = arr[midlatpx-int(lat_px_deg*8.5):midlatpx+int(lat_px_deg*8.5),midlongpx+int(long_px_deg*38):midlongpx+int(long_px_deg*57)]

plt.figure()
plt.imshow(arr_stereob[25:81,:],origin='lower',extent=(35,57,-8.5,8.5),aspect='auto')
plt.xlabel('Helioecliptic Longitude')
plt.ylabel('Helioecliptic Latitude')
plt.title('STEREO B')
plt.colorbar()
plt.savefig('res_disk_B.png')

xb = np.linspace(38,57,68)
pl = 600*np.power(xb,-1.6)

mean_b = np.average(arr_stereob[25:81,:],0)
plt.figure()
plt.plot(xb,mean_b,xb,pl)
plt.plot(xb,np.subtract(mean_b,pl))
plt.axhline(linestyle='dashed',color='0.75')
plt.xlim([35,57])
plt.savefig('mean_b.png')


arr_stereoa = arr[midlatpx-int(lat_px_deg*8.5):midlatpx+int(lat_px_deg*8.5),midlongpx-int(long_px_deg*57):midlongpx-int(long_px_deg*38)]
plt.figure()
plt.imshow(arr_stereoa[25:81,:],origin='lower',extent=(-57,-35,-8.5,8.5),aspect='auto')
plt.xlabel('Helioecliptic Longitude')
plt.ylabel('Helioecliptic Latitude')
plt.title('STEREO A')
plt.colorbar()
plt.savefig('res_disk_A.png')

xa = np.linspace(-57,-38,68)

mean_a = np.average(arr_stereoa[25:81,:],0)
plt.figure()
plt.plot(xa,mean_a)
plt.savefig('mean_a.png')

os.chdir('/home/yc367/Desktop/Disk_STEREO_images')