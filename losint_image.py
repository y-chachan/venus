from subprocess import Popen,PIPE
from scipy.io import readsav
from os import remove
from os.path import isfile

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
