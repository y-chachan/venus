#to take a 2D model and add thickness to it artifically
#then to look at a slice of that image
#to be used with Andrew's 2D models of circumsolar dust rings

#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import os
from astropy.convolution import Gaussian1DKernel
from matplotlib import gridspec
import pickle

os.chdir('/home/yc367/Desktop/Disk_Model')

arr = np.loadtxt('diskimage_beta=0.075.dat')
sur_den = arr[:,2].reshape((500,500))
resonance = arr[:,4].reshape((500,500))
size = (-1.8,1.8,-1.8,1.8)

plt.figure()
plt.imshow(resonance,origin='lower',extent=size)
plt.colorbar()
plt.savefig('resonance.png')

rel_angle = np.arange(0,361,15)   #the ange between Venus, Sun, and the tangent point 
#the ring moves with Venus around the Sun and therefore we need to 
#know how far away relative to Venus our tangent point is

img_height = 500     #image height in pixels
g = Gaussian1DKernel(stddev=100,x_size=500)

fig = plt.figure(figsize=(6,2*len(rel_angle)))
gs = gridspec.GridSpec(len(rel_angle),3)

for i in range(0,len(rel_angle)):
    rotated = ndimage.interpolation.rotate(resonance,rel_angle[i],(0,1))
    m_idx = int(len(rotated)/2)
    rot = rotated[(m_idx-250):(m_idx+250),(m_idx-250):(m_idx+250)]    
    rot_tD_image = np.zeros((img_height,np.shape(rot)[0],np.shape(rot)[1]))
    for j in range(0,img_height):
        rot_tD_image[j] = rot*(g.array[j])
    twoD_proj_x = np.sum(rot_tD_image,axis=1)
    twoD_proj_y = np.sum(rot_tD_image,axis=2)
    ax1 = fig.add_subplot(gs[i,0])
    ax1.imshow(rot,extent=size,origin='lower')
    ax1.set_ylabel('A.U.')
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = fig.add_subplot(gs[i,1])
    ax2.imshow(twoD_proj_x,extent=size,origin='lower')
    plt.setp(ax2.get_yticklabels(), visible=False)    
    plt.setp(ax2.get_xticklabels(), visible=False)    
    ax3 = fig.add_subplot(gs[i,2])
    ax3.imshow(twoD_proj_y,extent=size,origin='lower')
    plt.setp(ax3.get_xticklabels(), visible=False)    
    plt.setp(ax3.get_yticklabels(), visible=False)

plt.setp(ax1.get_xticklabels(),visible=True)
plt.setp(ax2.get_xticklabels(),visible=True)
plt.setp(ax3.get_xticklabels(),visible=True)
plt.setp([a.locator_params(nbins=6) for a in fig.axes])
ax1.set_xlabel('A.U.')
ax2.set_xlabel('A.U.')
ax3.set_xlabel('A.U.')
plt.savefig('disk_view.png')

for i in range(0,len(rel_angle)):
    rotated = ndimage.interpolation.rotate(resonance,rel_angle[i],(0,1))
    m_idx = int(len(rotated)/2)
    rot = rotated[(m_idx-250):(m_idx+250),(m_idx-250):(m_idx+250)]    
    rot_tD_image = np.zeros((img_height,np.shape(rot)[0],np.shape(rot)[1]))
    for j in range(0,img_height):
        rot_tD_image[j] = rot*(g.array[j])
    twoD_proj_x = np.sum(rot_tD_image,axis=1)
    twoD_proj_y = np.sum(rot_tD_image,axis=2)
    fig = plt.figure(figsize=(6,2))
    gs = gridspec.GridSpec(1,3)
    ax1 = fig.add_subplot(gs[0,0])
    ax1.imshow(rot,extent=size,origin='lower')
    ax1.set_ylabel('A.U.')
    ax1.set_xlabel('A.U.')
    ax2 = fig.add_subplot(gs[0,1])
    ax2.imshow(twoD_proj_x,extent=size,origin='lower')
    plt.setp(ax2.get_yticklabels(), visible=False)    
    ax2.set_xlabel('A.U.')
    ax3 = fig.add_subplot(gs[0,2])
    ax3.imshow(twoD_proj_y,extent=size,origin='lower')
    ax3.set_xlabel('A.U.')   
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp([a.locator_params(nbins=6) for a in fig.axes])
    plt.savefig('disk_view'+str(i)+'.png')
    plt.close()     #to prevent all the images from being displayed"""