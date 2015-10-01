from pylab import*
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import PyGuide
from scipy import ndimage#from scipy.interpolate import shift
import csv
import glob

#help(PyGuide.centroid)
path_12="/Users/jessicaluna/Documents/UTAustin/PlanetaryAstrophysics/Homework2/KOA_28560/NIRC2/calibrated/"
path_42B="/Users/jessicaluna/Documents/UTAustin/PlanetaryAstrophysics/Homework2/KOA_29864/NIRC2/calibrated"


#help(ndimage.interpolation.shift)


files_12=glob.glob(path_12 + '/*.fits')
files_42B=glob.glob(path_42B+'/*.fits')

len_12=len(files_12)
len_42B=len(files_42B)

testfile=fits.open(files_42B[0])
#help(test)
testfile_data=testfile[0].data

#### part 2 finding x and y pixels of flux weighted centroid


centroid_xlist_12=[]
centroid_ylist_12=[]

centroid_xlist_42B=[]
centroid_ylist_42B=[]

PA_42B=[]
PA_12=[]
##for ROXs 42B
for indexb in range (0,len_42B):
    imageb=fits.open(files_42B[indexb])
    #header=imagetest[0].header['filename']
    countposxb=[]
    countposyb=[]
    datab=imageb[0].data
    datab1=datab[0]
    datab2=datab[1]
    for iii in range (605,620):
        countposb=datab1[iii]*iii
        countposxb.append(countposb)
    for jjj in range (462,477):
        countposc=datab2[jjj]*jjj
        countposyb.append(countposc)
    totalcountsxb=sum(datab1[605:620])
    totalcountsyb=sum(datab2[462:477])
    sumcountposxb=sum(countposxb)
    sumcountposyb=sum(countposyb)
    centroid_xb=sumcountposxb/totalcountsxb
    centroid_yb=sumcountposyb/totalcountsyb
    centroid_xlist_42B.append(centroid_xb)
    centroid_ylist_42B.append(centroid_yb)

    PARANG_a= imageb[0].header['PARANG']
    ROTPPOSN_a= imageb[0].header['ROTPPOSN']
    EL_a= imageb[0].header['EL']
    INSTANGL_a= imageb[0].header['INSTANGL']
    
    PA_a= PARANG_a + ROTPPOSN_a - EL_a - INSTANGL_a
    PA_42B.append(PA_a)

## for ROXS 12
## the coronagraph seems to be centered on a box from 600-620 in x and 460-480 in y
for index in range (0,len_12):
    imagetest=fits.open(files_12[index])
    #header=imagetest[0].header['filename']
    countposx=[]
    countposy=[]
    data=imagetest[0].data
    data1=data[0]
    data2=data[1]
    for ii in range (605,620):
        countposa=data1[ii]*ii
        countposx.append(countposa)
    for jj in range (462,477):
        countposb=data2[jj]*jj
        countposy.append(countposb)
    totalcountsx=sum(data1[605:620])
    totalcountsy=sum(data2[462:477])
    sumcountposx=sum(countposx)
    sumcountposy=sum(countposy)
    centroid_x=sumcountposx/totalcountsx
    centroid_y=sumcountposy/totalcountsy
    centroid_xlist_12.append(centroid_x)
    centroid_ylist_12.append(centroid_y)


    PARANG_b= imagetest[0].header['PARANG']
    ROTPPOSN_b= imagetest[0].header['ROTPPOSN']
    EL_b= imagetest[0].header['EL']
    INSTANGL_b= imagetest[0].header['INSTANGL']
    
    PA_b= PARANG_b + ROTPPOSN_b - EL_b - INSTANGL_b
    PA_12.append(PA_b)
    
savetxt('centroids_roxs12.txt', zip(centroid_xlist_12,centroid_ylist_12), fmt='%.18e', delimiter=' ', newline='\n', header='centroids ROXs12')

savetxt('centroids_roxs42b.txt', zip(centroid_xlist_42B,centroid_ylist_42B), fmt='%.18e', delimiter=' ', newline='\n', header='centroids ROXs42b')

savetxt('PA_42B.txt',PA_42B, fmt='%.18e', delimiter=' ', newline='\n', header='PAS 42B')

savetxt('PA_12.txt',PA_12, fmt='%.18e', delimiter=' ', newline='\n', header='PAS 12')

