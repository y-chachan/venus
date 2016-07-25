from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import optimize
from matplotlib import gridspec
import pickle
import functions
import glob

path1 = 'C:/Users/YaChan/Desktop/Output_A'   #output folder
path = 'Y:/IoA Summer Project/STEREO/A'
fold_names = [x[0] for x in os.walk(path)]
length = len(fold_names)
l1 = 6
for l1 in range(6,7):    
    path2 = fold_names[l1]   #input data
    output_file = os.path.relpath(path2,start='Y:/IoA Summer Project/STEREO/A')
    os.chdir(path2)
    
    #import all the fits data from the folder
    fnames = glob.glob('*.fts')
    file_count = len(fnames)
    
    #header information
    header = fits.getheader(fnames[0]) #print(header['CDELT1']) #example
    
    pxsize = header['CDELT1']   #angular size of a pixel
    midln = header['CRVAL1']    #the longitude at the centre of the image
    midlt = header['CRVAL2']    #the latitude at the centre of the image
    size = header['NAXIS1']     #number of pixels, we only deal with 1024x1024 images
    var = size/2                
    
    #selecting the relevant area of the images (the numerals added at the end
    #are to make sure we can divide the pixels into an integer number of cells)
    ilng = -55  #to change the HE longitude, also DIFFERENT for A (-60) and B (35)
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
    
    j = 0
    for file in fnames:
        hdulist = fits.open(file)
        dc[j] = hdulist[0].data
        header = hdulist[0].header
        hdulist.close()
        time[j] = header['DATE-OBS']
        stereo_lon[j] = np.degrees(np.arctan(header['HAEY_OBS']/header['HAEX_OBS']))
        data_cube[j] = dc[j,lt1:lt2,ln1:ln2]
        j += 1
    del j, var
    
    #reshape and create an 2d array from all the different images
    rs_data_cube = np.swapaxes(data_cube,0,2).reshape((pix_deg*hcellsize*hgridsize,vgridsize,pix_deg*vcellsize*file_count)).swapaxes(0,1).reshape(vgridsize,hgridsize,pix_deg*hcellsize*pix_deg*vcellsize*file_count)
    image = np.zeros((vgridsize,hgridsize)) #to be used to store histogram values
    
    #calculate a certain medial value and assign the entire cell that value
    medianvalue = 0.45      #0.45 value (55th percentile)
    for a in range(0,hgridsize):
        for b in range(0,vgridsize):
            image[b,a] = np.nanpercentile(rs_data_cube[b,a,:],medianvalue*100)
    del a,b
         
    # define your function for leastsq to obtain initial value estimates:
    def func(x):
        return k() * np.power(np.absolute(x),n())
      
    IV = np.zeros((vgridsize,2))
    k = functions.Parameter(10)
    n = functions.Parameter(-1.5)
    for i in range(0,vgridsize):
        IV[i] = functions.fit(func,[k,n],image[i,:],lg)[0]   #gives the intial values we can use for the more accurate method
    del i 
    # define your function: for fitting using curve_fit
    def f(x,k,n,c): 
        return (k*np.power(np.absolute(x),n)) + c 
        
    # providing the data, fitting power laws for constant beta (he-latitude)
    number_parameters = 3
    c = np.zeros(vgridsize)
    
    #extrapolate power law data
    #pl_image = np.zeros((vgridsize,hgridsize))
    
    for i in range(0,vgridsize):
        k = functions.Parameter(IV[i,0])
        n = functions.Parameter(IV[i,1])
        #pl_image[i] = func(lg)
        c[i] = func(lg[1])-image[i,1]
    del i
    power_law = np.zeros((vgridsize,number_parameters))
    for i in range(0,vgridsize):
        power_law[i] = optimize.curve_fit(f,lg,image[i,:],[IV[i,0],IV[i,1],c[i]])[0]
    del i
    #extrapolate power law data
    pl_image = np.zeros((vgridsize,hgridsize))
    
    for i in range(0,vgridsize):
        pl_image[i] = f(lg,power_law[i,0],power_law[i,1],power_law[i,2])
    del i
      
    #detrending data and power law extrapolation using boxcar filtering
    filter_width = 13       #assign the filter width in pixels
    data_filtered = np.zeros((vgridsize,hgridsize))
    for i in range(0,vgridsize):
        data_filtered[i] = convolve(image[i,:],Box1DKernel(filter_width))
    del i
    smooth_pl = np.zeros((vgridsize,hgridsize))
    for i in range(0,vgridsize):
        smooth_pl[i] = convolve(pl_image[i,:],Box1DKernel(filter_width)) 
    del i
    
    os.chdir(path1)     #changing directory to Output before we start saving figures
    
    #subtracting the two results and obtaining the final result
    inter_image = np.subtract(image,data_filtered)
    inter_pl = np.subtract(pl_image,smooth_pl)
    
    final_image = np.subtract(inter_image,inter_pl)
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
        ax3.plot(lg,np.add(final_image[i,:],i/5))
        ax3.hlines(i/5,ilng,flng,linestyles='dashed',color = '0.75')
    ax3.set_ylim([-0.1,3.5])
    ax3.set_xlim([ilng,flng])
    ax3.set_title('Residual')
    ax3.set_xlabel('Helioecliptic Longitude')
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0)
    plt.savefig(output_file+'.png')
    plt.close()
    del i
    
    fig3 = plt.figure(figsize=(8,vgridsize*2))
    gs3 = gridspec.GridSpec(vgridsize,2)
    for i in range(0,vgridsize):
        pl = np.power(np.absolute(lg),abs(IV[i,1]))    
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
    plt.setp(ax3.get_xticklabels(),visible=True)
    #gs3.tight_layout(fig3)
    plt.savefig(output_file+'_diagnostics.png')
    plt.close()
    del i
    l1 += 1