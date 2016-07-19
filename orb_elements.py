from __future__ import division 	#necessary in older versions of Python
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec

stereo_a = np.genfromtxt('stereoa_orb.txt',skip_header=24,skip_footer=40,delimiter=",")
stereo_b = np.genfromtxt('stereob_orb.txt',skip_header=24,skip_footer=40,delimiter=",")
time = np.linspace(2008,2016+(9/12),98,endpoint=True)

a_e = np.divide(np.subtract(stereo_a[:,2],np.mean(stereo_a[:,2])),stereo_a[:,2])
a_inc = np.divide(np.subtract(stereo_a[:,4],np.mean(stereo_a[:,4])),stereo_a[:,4])
a_asc = np.divide(np.subtract(stereo_a[:,5],np.mean(stereo_a[:,5])),stereo_a[:,5])
a_per = np.divide(np.subtract(stereo_a[:,6],np.mean(stereo_a[:,6])),stereo_a[:,6])
a_a = np.divide(np.subtract(stereo_a[:,11],np.mean(stereo_a[:,11])),stereo_a[:,11])

b_e = np.divide(np.subtract(stereo_b[:,2],np.mean(stereo_b[:,2])),stereo_b[:,2])
b_inc = np.divide(np.subtract(stereo_b[:,4],np.mean(stereo_b[:,4])),stereo_b[:,4])
b_asc = np.divide(np.subtract(stereo_b[:,5],np.mean(stereo_b[:,5])),stereo_b[:,5])
b_per = np.divide(np.subtract(stereo_b[:,6],np.mean(stereo_b[:,6])),stereo_b[:,6])
b_a = np.divide(np.subtract(stereo_b[:,11],np.mean(stereo_b[:,11])),stereo_b[:,11])

fig = plt.figure(figsize=(8,10))
gs = gridspec.GridSpec(2,1)
ax1 = fig.add_subplot(gs[0,0])
e, = ax1.plot(time,a_e)
inc, = ax1.plot(time,a_inc)
asc, = ax1.plot(time,a_asc)
per, = ax1.plot(time,a_per)
a, = ax1.plot(time,a_a)
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
ax2 = fig.add_subplot(gs[1,0])
ax2.plot(time,b_e)
ax2.plot(time,b_inc)
ax2.plot(time,b_asc)
ax2.plot(time,b_per)
ax2.plot(time,b_a)
ax1.set_ylim([-0.02,0.02])
ax2.set_ylim([-0.005,0.005])
ax1.set_ylabel('Fractional change in Orbital Elements')
ax1.set_xlabel('Time')
ax2.set_ylabel('Fractional change in Orbital Elements')
ax2.set_xlabel('Time')
ax1.set_title('STEREO A')
ax2.set_title('STEREO B')
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.legend([e,inc,asc,per,a], ['Eccentricity', 'Inclination', 'Ascending Node','Pericenter Long.','Semi-major Axis'],bbox_to_anchor=(0, 2.35, 1., .102), loc=9,ncol=3, borderaxespad=0.)
plt.savefig('orb_elements.png')