from astropy.io import fits
import numpy as np
import os
import pickle
from functions import *
import glob

path1 = 'C:/Users/YaChan/Desktop'   #output folder
path2 = 'Y:/IoA Summer Project/STEREO/B'    #input folder
fold_names = [x[0] for x in os.walk(path)]
length = len(fold_names)
l1 = 2
for l1 in range(2,3):    
    path = fold_names[l1]   #input data
    output_file = os.path.relpath(path,start=path2)
    os.chdir(path2)
    
    #import all the fits data from the folder
    fnames = glob.glob('*.fts')
    file_count = len(fnames)
    
    #header information
    header = fits.getheader(fnames[0]) #print(header['CDELT1']) #example
    sat = header['OBSRVTRY']
    pxsize = header['CDELT1']   #angular size of a pixel
    midln = header['CRVAL1']    #the longitude at the centre of the image
    midlt = header['CRVAL2']    #the latitude at the centre of the image
    size = header['NAXIS1']     #number of pixels, we only deal with 1024x1024 images
    var = size/2                
    
    #selecting the relevant area of the images (the numerals added at the end
    #are to make sure we can divide the pixels into an integer number of cells)
    if (sat == 'STEREO_A')  :
        ilng = -60
        flng = -35
    elif (sat == 'STEREO_B')  :
        ilng = 35
        flng = 60
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
    
    output = process_image(image,vgridsize,hgridsize,lg)    
    
    os.chdir(path1)     #changing directory to Output before we start saving figures
    
    #subtracting the two results and obtaining the final result
    inter_image = np.subtract(image,output[0])
    inter_pl = np.subtract(output[1],output[2])
    final_image = np.subtract(inter_image,inter_pl)
    
    #saving important output to a pickle file
    pkl_file = output_file+'.p'
    data = {
    'image': image,
    'pl' : output[3],
    'final_image' : final_image
    }    
    pickle.dump(data,open(pkl_file, "wb" ))
    
    plot_output(final_image,vgridsize,lg,ilng,flng,ilat,flat,output_file)
    plot_diagnostics(image,output[0],output[1],output[2],inter_image,inter_pl,final_image,vgridsize,lg,ilng,flng,output[3],output_file)
    
    l1 += 1