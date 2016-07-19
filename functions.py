import numpy as np
from scipy import optimize
import os, glob
from astropy.io import fits
import matplotlib.pyplot as plt

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
    