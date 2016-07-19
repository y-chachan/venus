from astropy.io import fits
#from astropy.convolution import convolve, Box1DKernel
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import optimize
from matplotlib import gridspec
#import operator
import pickle
#from astropy.convolution import Gaussian1DKernel
import functions
import glob

path1 = 'C:/Users/YaChan/Desktop/Output'   #output folder
path2 = 'Y:/IoA Summer Project/STEREO/A/20131201'    #input data
output_file = os.path.relpath(path2,start='Y:/IoA Summer Project/STEREO/A')
os.chdir(path2)

#import all the fits data from the folder
fnames = glob.glob('*.fts')
file_count = len(fnames)

#header information
header = fits.getheader(fnames[0])  #print(header['CDELT1']) #example

pxsize = header['CDELT1']   #angular size of a pixel
midln = header['CRVAL1']    #the longitude at the centre of the image
midlt = header['CRVAL2']    #the latitude at the centre of the image
size = header['NAXIS1']     #number of pixels, we only deal with 1024x1024 images
var = size/2                

#selecting the relevant area of the images (the numerals added at the end
#are to make sure we can divide the pixels into an integer number of cells)
ilng = -60  #to change the HE longitude, also DIFFERENT for A (-60) and B (35)
flng = -35  #to change the HE longitude, also DIFFERENT for A (-35) and B (60)
ilat = -8.5
flat = 8.5
hextent = flng-ilng    #horizontal extent of the image (35 to 60 degrees)
vextent = flat-ilat    #vertical extent(-8.5 to 8.5 degrees)
hcellsize = 0.5        #cell size in the horizontal direction (degrees)
vcellsize = 1          #cell size in the vertical direction (degrees)
pix_deg = round(1/pxsize)    #cell size (1 degree in pixels)
hgridsize = int(hextent/hcellsize)   #size of the recast images
vgridsize = int(vextent/vcellsize)   #size of the recast images

lg = (hcellsize)*np.arange(hgridsize) + ilng
lt = (vcellsize)*np.arange(vgridsize) + ilat

lt1 = var - round((midlt - ilat)*pix_deg)
lt2 = var + round((flat - midlt)*pix_deg)
ln1 = var - round((midln - ilng)*pix_deg)
ln2 = var + round((flng - midln)*pix_deg)

dc = np.zeros((file_count,size,size))
data_cube = np.zeros((file_count,lt2-lt1,ln2-ln1))

time = np.zeros(file_count).astype(object)
stereo_lon = np.zeros(file_count)

i = 0
for file in fnames:
    hdulist = fits.open(file)
    dc[i] = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    time[i] = header['DATE-OBS']
    stereo_lon[i] = np.degrees(np.arctan(header['HAEY_OBS']/header['HAEX_OBS']))
    data_cube[i] = dc[i,lt1:lt2,ln1:ln2]
    i += 1
del i, var

#reshape and create an 2d array from all the different images
rs_data_cube = np.swapaxes(data_cube,0,2).reshape((pix_deg*hcellsize*hgridsize,vgridsize,pix_deg*vcellsize*file_count)).swapaxes(0,1).reshape(vgridsize,hgridsize,pix_deg*hcellsize*pix_deg*vcellsize*file_count)
image = np.zeros((vgridsize,hgridsize)) #to be used to store histogram values

#calculate a certain medial value and assign the entire cell that value
medianvalue = 0.45      #0.45 value (55th percentile)
for j in range(0,hgridsize):
    for i in range(0,vgridsize):
        image[i,j] = np.nanpercentile(rs_data_cube[i,j,:],medianvalue*100)
del i,j

# define your function for leastsq to obtain initial value estimates:
def func(x):
    return k() * np.power(np.absolute(x),n())
  
IV = np.zeros((vgridsize,2))
k = functions.Parameter(1)
n = functions.Parameter(-1.5)
for i in range(0,vgridsize):
    IV[i] = functions.fit(func,[k,n],image[i,:],lg)[0]   #gives the intial values we can use for the more accurate method

# define your function: for fitting using curve_fit
def f(x,k,n,c): 
    return (k*np.power(np.absolute(x),n)) + c 
    
# providing the data, fitting power laws for constant beta (he-latitude)
number_parameters = 3
c = np.zeros(vgridsize)
for i in range(0,vgridsize):
    k = functions.Parameter(IV[i,0])
    n = functions.Parameter(IV[i,1])
    c[i] = func(lg[0])-image[i,0]

power_law = np.zeros((vgridsize,number_parameters))
for i in range(0,vgridsize):
    power_law[i] = optimize.curve_fit(f,lg,image[i,:],[IV[i,0],IV[i,1],c[i]])[0]

#extrapolate power law data
pl_image = np.zeros((vgridsize,hgridsize))

for i in range(0,vgridsize):
    pl_image[i] = f(lg,power_law[i,0],power_law[i,1],power_law[i,2])

"""    
#detrending data and power law extrapolation using boxcar filtering
filter_width = 1       #assign the filter width in pixels
data_filtered = np.zeros((vgridsize,hgridsize))
for i in range(0,vgridsize):
    data_filtered[i] = convolve(image[i,:],Box1DKernel(filter_width))

smooth_pl = np.zeros((vgridsize,hgridsize))
for i in range(0,vgridsize):
    smooth_pl[i] = convolve(pl_image[i,:],Box1DKernel(filter_width)) 
del i
"""
os.chdir(path1)     #changing directory to Output before we start saving figures

#subtracting the two results and obtaining the final result
final_image = np.subtract(image,pl_image)
#mean scan (add all values for a given longitude between -4.5 to 4.5 degrees latitude)
mean_scan = np.average(final_image,0)

pkl_file = output_file+'.p'
pickle.dump(final_image,open(pkl_file, "wb" ))

fig = plt.figure(figsize=(7,6))
gs = gridspec.GridSpec(5,6)
ax1 = fig.add_subplot(gs[0,0:4])
ax1.plot(lg[6:44],mean_scan[6:44])
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
    ax3.plot(lg,np.add(np.subtract(image[i,:],pl_image[i,:]),i/5))
    ax3.hlines(i/5,ilng,flng,linestyles='dashed',color = '0.75')
ax3.set_ylim([-0.1,3.5])
ax3.set_xlim([ilng,flng])
ax3.set_title('Residual')
ax3.set_xlabel('Helioecliptic Longitude')
plt.tight_layout()
#plt.subplots_adjust(hspace=0)
plt.savefig(output_file+'.png')