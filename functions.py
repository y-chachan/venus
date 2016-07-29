from __future__ import division
import numpy as np
from scipy import optimize
import os, glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from matplotlib import gridspec
import ephem
from subprocess import Popen,PIPE
from scipy.io import readsav
from os import remove
from os.path import isfile

def process_image(medianed_image,vgridsize,hgridsize,lg,filter_width=13):
    # define your function for leastsq to obtain initial value estimates:
    def func(x):
        return k() + np.multiply(n(),np.absolute(x))
        
    IV = np.zeros((vgridsize,2))
    k = Parameter(5)
    n = Parameter(-1.5)
    for i in range(0,vgridsize):
        IV[i] = fit(func,[k,n],np.log10(medianed_image[i,:]),np.log10(np.abs(lg)))[0]

    #extrapolate power law data
    pl_image = np.zeros((vgridsize,hgridsize))    
    for i in range(0,vgridsize):
        k = Parameter(IV[i,0])
        n = Parameter(IV[i,1])
        pl_image[i] = (10**IV[i,0])*np.power(np.abs(lg),IV[i,1])    
      
    #detrending data and power law extrapolation using boxcar filtering
    data_filtered = np.zeros((vgridsize,hgridsize))
    for i in range(0,vgridsize):
        data_filtered[i] = convolve(medianed_image[i,:],Box1DKernel(filter_width))
    
    smooth_pl = np.zeros((vgridsize,hgridsize))
    for i in range(0,vgridsize):
        smooth_pl[i] = convolve(pl_image[i,:],Box1DKernel(filter_width)) 
    del i
    return data_filtered,pl_image,smooth_pl,IV

def plot_output(final_image,vgridsize,lg,ilng,flng,ilat,flat,output_file,filter_width=13):
    mean_scan = np.average(final_image,0)
    cut = int(filter_width/2)
    fig = plt.figure(figsize=(7,6))
    gs = gridspec.GridSpec(5,6)
    ax1 = fig.add_subplot(gs[0,0:4])
    ax1.plot(lg[cut:-cut],mean_scan[cut:-cut])
    ax2 = fig.add_subplot(gs[1:5,0:4],sharex=ax1)
    ax1.set_title('Mean Scan')
    ax1.set_ylabel('DN s$^{-1}$pix$^{-1}$')
    ax1.hlines(0,ilng,flng,linestyles='dashed',color = '0.75')
    plt.setp(ax1.get_xticklabels(), visible=False)
    im = ax2.imshow(final_image[:,6:44],origin='lower',aspect='auto',extent=(ilng+3,flng-3,ilat,flat))
    ax2.set_xlabel('Helioecliptic Longitude')
    ax2.set_ylabel('Helioecliptic Latitude')
    plt.colorbar(im,ax=ax2,orientation='horizontal')
    ax3 = fig.add_subplot(gs[0:5,4:6])
    for i in range(0,vgridsize):
        ax3.plot(lg,np.add(final_image[i,:],i/5))
        ax3.hlines(i/5,ilng,flng,linestyles='dashed',color = '0.75')
    ax3.set_ylim([-0.1,vgridsize*0.2])
    ax3.set_xlim([ilng,flng])
    ax3.set_title('Residual')
    ax3.set_xlabel('Helioecliptic Longitude')
    plt.tight_layout()
    plt.savefig(output_file+'.png')
    plt.close()

def plot_diagnostics(image,data_filtered,pl_image,smooth_pl,inter_image,inter_pl,final_image,vgridsize,lg,ilng,flng,IV,output_file):
    fig3 = plt.figure(figsize=(8,vgridsize*2))
    gs3 = gridspec.GridSpec(vgridsize,2)
    for i in range(0,vgridsize):
        pl = np.power(lg,abs(IV[i,1]))    
        ax1 = fig3.add_subplot(gs3[i,0])
        [a,b,c,d] = ax1.plot(lg,np.multiply(pl,image[i,:]),lg,np.multiply(pl,data_filtered[i,:]),lg,np.multiply(pl,pl_image[i,:]),lg,np.multiply(pl,smooth_pl[i,:]))
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_ylim([min(np.multiply(pl,image[i,:])),max(np.multiply(pl,image[i,:]))]) 
        plt.locator_params(axis='y',nbins=5)    
        if i < 1:
            ax1.set_title('Raw Images')
            plt.legend([a,b,c,d], ['Data', 'Data boxcar', 'Power Law','Power Law boxcar'],bbox_to_anchor=(0, 2.0, 1., .102),loc=9,
               borderaxespad=0.)
        ax2 = fig3.add_subplot(gs3[i,1],sharex=ax1)
        [e,f,g] = ax2.plot(lg,inter_image[i,:],lg,inter_pl[i,:],lg,final_image[i,:])
        plt.setp(ax2.get_xticklabels(), visible=False) 
        ax2.set_xlim([ilng+3,flng-3])
        ax2.set_ylim([min(inter_image[i,:]),max(final_image[i,:])])     
        ax2.hlines(0,ilng,flng,linestyles='dashed',color='0.75')
        #plt.locator_params(axis='y',nbins=5)     
        if i <1 :
            ax2.set_title('Subtracted Images')
            plt.legend([e,f,g], ['Subtracted Image', 'Subtracted Power Law', 'Final Image'],bbox_to_anchor=(0, 1.8, 1., .102),loc=9,
               borderaxespad=0.)
    plt.setp(ax1.get_xticklabels(),visible=True)
    plt.setp(ax2.get_xticklabels(),visible=True)
    plt.savefig(output_file+'_diagnostics.png')
    plt.close()

def diagnose(path1,path2):
    """
    To check if given image(s) can be used to detect the circumsolar ring
    Args:
      path1: the directory from which a reference image is created
      path2: the collection of images that need to be diagnosed
    """
    #create a reference image
    #path1 = 'C:/Users/YaChan/Desktop/STEREO/B/20090629'
    os.chdir(path1)
    fnames = glob.glob('*.fts')
    size = 1024
    cell = 16
    grid = int(size/cell)    
    file_count = len(fnames)
    dc = np.zeros((file_count,size,size))
    for j in range(0,file_count):
        dc[j] = fits.getdata(fnames[j])
    rs_data_cube = np.swapaxes(dc,0,2).reshape((size,grid,file_count*cell)).swapaxes(0,1).reshape(grid,grid,cell*cell*file_count)
    ref_image = np.zeros((grid,grid))
    medianvalue = 0.45
    for k in range(0,grid):
        for l in range(0,grid):
                ref_image[k,l] = np.nanpercentile(rs_data_cube[k,l,:],medianvalue*100)
    
    ref_median = np.median(ref_image)

    #diagnose the supplied image
    os.chdir(path2)
    file_names = glob.glob('*.fts')
    #header = fits.getheader(file_names[0])   
    file_count = len(file_names)
    dc = np.zeros((file_count,size,size))
    for j in range(0,file_count):
        dc[j] = fits.getdata(file_names[j])
    rs_data_cube = np.swapaxes(dc,0,2).reshape((size,grid,file_count*cell)).swapaxes(0,1).reshape(grid,grid,cell*cell*file_count)
    image = np.zeros((grid,grid))
    medianvalue = 0.45
    for k in range(0,grid):
        for l in range(0,grid):
            image[k,l] = np.nanpercentile(rs_data_cube[k,l,:],medianvalue*100)
    residue_median = np.median(image) - ref_median
    plt.figure()
    plt.imshow(np.subtract(image,ref_image),origin='lower')
    plt.title('Given Image subtracted by Reference Image')    
    plt.colorbar()
    print('Residue median =' + residue_median)


class Parameter:
    def __init__(self, value):
            self.value = value

    def set(self, value):
            self.value = value

    def __call__(self):
            return self.value

def fit(function, parameters, y, x):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)

    p = [param() for param in parameters]
    return optimize.leastsq(f, p, maxfev=1000)
    

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
