from pylab import*
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import ndimage#from scipy.interpolate import shift
import csv
import glob
import time
#from photutils import daofind

starttime=time.time()

path_12="/Users/jessicaluna/Documents/UTAustin/PlanetaryAstrophysics/Homework2/KOA_28560/NIRC2/calibrated/"
path_42B="/Users/jessicaluna/Documents/UTAustin/PlanetaryAstrophysics/Homework2/KOA_29864/NIRC2/calibrated"




files_12=glob.glob(path_12 + '/*.fits')
files_42B=glob.glob(path_42B+'/*.fits')

len_12=len(files_12)
len_42B=len(files_42B)


#### part 2 importing centroid lists from centroid.py



centroid_xlist_42B,centroid_ylist_42B   =loadtxt('centroids_roxs42b.txt', unpack='True', usecols=[0,1])
centroid_xlist_12,centroid_ylist_12   =loadtxt('centroids_roxs12.txt', unpack='True', usecols=[0,1])

PA_42B=loadtxt('PA_42b.txt', unpack='True', usecols=[0])
PA_12=loadtxt('PA_12.txt', unpack='True', usecols=[0])


newdata=[]
rotated_42B=[]
subtract_median=[]


for ind in range (0,len_42B):
    #shifting all star position to center of images
    shiftx= 512.0 - centroid_xlist_42B[ind]
    shifty= 512.0 - centroid_ylist_42B[ind]
            
    image1=fits.open(files_42B[ind])
    data_image1=image1[0].data
    shifts_42 = ndimage.interpolation.shift(data_image1,[shifty ,shiftx])
    shifts_42_profile= ndimage.interpolation.shift(data_image1,[shifty ,shiftx])
    
    newdata.append(shifts_42)

    #rotation about star(center of image)

    imgR = ndimage.rotate(shifts_42, -PA_42B[ind] , reshape=False)

    rotated_42B.append(imgR)
    
    #part 5:
    data_subtract= shifts_42_profile
    yy=linspace(0,1024,1024)
    xx=linspace(0,1024,1024)
    y,x=np.meshgrid(yy,xx)
    xcord= x- 512.0
    ycord= y - 512.0
    radius= sqrt(xcord**2+ycord**2)
    #radii.append(radius)
    
    #print(len(xcord))
    #print(len(radius))
    ring_len=2.0
    for rr in range(0,100):
        ind1=np.where( (radius>rr*ring_len) & (radius<(rr+1)*ring_len) )
        ring1=data_subtract[ind1]
        median1=np.median(ring1)
        data_subtract[ind1]=data_subtract[ind1]-median1
    subtract_median.append(data_subtract)

Median_images_42B=np.median(newdata,axis=0)

sum_42B=np.sum(newdata,axis=0)


rotated_median_42B=np.median(rotated_42B, axis=0)
rotated_sum_42B=np.sum(rotated_42B,axis=0)

med_rot_subtract=np.median(subtract_median,axis=0)
sum_rot_subtract=np.sum(subtract_median,axis=0)
           




#part 6 subtract median PSF from each image
ADI_images_42B=[]
for indi in range (0,len_42B):
    ADI_imagesa= newdata[indi] - Median_images_42B
    ADI_images_42B.append(ADI_imagesa)

ADI_med_42B=np.median(ADI_images_42B,axis=0)
ADI_sum_42B=np.sum(ADI_images_42B,axis=0)



        
newdata_12=[]
rotated_12=[]
subtract_median_12=[]
for ind_12 in range (0,len_12):
    ##shifting all images to center of image
    shiftx_12=512.0-centroid_xlist_12[ind_12]
    shifty_12= 512.0-centroid_ylist_12[ind_12]
    
        
    image1_12=fits.open(files_12[ind_12])
    data_image1_12=image1_12[0].data
    
    shifts_12 = ndimage.interpolation.shift(data_image1_12,[shifty_12 ,shiftx_12])
    shifts_12_profile = ndimage.interpolation.shift(data_image1_12,[shifty_12 ,shiftx_12])
    
    newdata_12.append(shifts_12) 


    
    imgRb = ndimage.rotate(shifts_12, -PA_12[ind_12] , reshape=False)
    rotated_12.append(imgRb)

 part 5
    data_subtract_12= shifts_12_profile
    yya=linspace(0,1024,1024)
    xxa=linspace(0,1024,1024)
    ya,xa=np.meshgrid(yya,xxa)
    xcorda= xa- 512.0
    ycorda= ya - 512.0
    radiusa= sqrt(xcorda**2+ycorda**2)
    #radii.append(radius)
    
    #print(len(xcord))
    #print(len(radius))
    ring_len=2.0
    for rrr in range(0,100):
        ind1a=np.where( (radiusa>rrr*ring_len) & (radiusa<(rrr+1)*ring_len) )
        ring1a=data_subtract_12[ind1a]
        median1a=np.median(ring1a)
        data_subtract_12[ind1a]=data_subtract_12[ind1a]-median1a
    subtract_median_12.append(data_subtract_12)

    
sum_12=np.sum(newdata_12,axis=0)
Median_images_12=np.median(newdata_12,axis=0)

rotated_median_12=np.median(rotated_12, axis=0)
rotated_sum_12=np.sum(rotated_12,axis=0)

med_rot_subtract_12=np.median(subtract_median_12,axis=0)
sum_rot_subtract_12=np.sum(subtract_median_12,axis=0)


#part 6 subtract median PSF from each image
ADI_images_12=[]
for indj in range (0,len_12):
    ADI_imagesb= newdata_12[indj] - Median_images_12
    ADI_images_12.append(ADI_imagesb)

ADI_med_12=np.median(ADI_images_12, axis=0)
ADI_sum_12=np.sum(ADI_images_12, axis=0)


#part 7 comparing each image to find the one with the most similar PSF
part_7=[]
file_most_sim=[]
for uu in range (0, len_42B):
    diff_42B=[]
    img=newdata[uu]
    
    for nn in range(0,len_12):
        diff=img-newdata_12[nn]
        difftot=np.sum(abs(diff))
    #print(difftot)
        diff_42B.append(difftot)
#diff_42Ba=np.asarray(diff_42B)  
    min_index=diff_42B.index(min(diff_42B))
#min_index=np.where( diff_42B == mina  )

#print(len(min_index))
    #print(min_index)
    file_most_sim.append(files_12[min_index])
    bestpsf= img-newdata_12[min_index]
    part_7.append(bestpsf)
    
#print(file_most_sim)

med_part7=np.median(part_7,axis=0)
sum_part7=np.sum(part_7,axis=0)

part_7a=[]
file_most_sima=[]
for uuu in range (0, len_12):
    diff_12=[]
    imga=newdata_12[uuu]
    
    for nnn in range(0,len_42B):
        diffa=imga-newdata[nnn]
        difftota=np.sum(abs(diffa))
    #print(difftot)
        diff_12.append(difftota)
#diff_42Ba=np.asarray(diff_42B)  
    min_indexa=diff_12.index(min(diff_12))
#min_index=np.where( diff_42B == mina  )

#print(len(min_index))
    #print(min_indexa)
    file_most_sima.append(files_42B[min_indexa])
    bestpsfa= imga-newdata[min_indexa]
    part_7a.append(bestpsfa)
    
#print(file_most_sima)
#print(file_most_sima[35], file_most_sim[6])
med_part7a=np.median(part_7a,axis=0)
#sum_part7a=np.sum(part_7a,axis=0)



#finding position of objects by eye
arcsec_pixel=0.009942

## center for part 7 roxs42B
xdist_742=abs(512.0-622.0)
ydist_742=abs(512.0-473.5)

##center for part 7 roxs12
xdist_712=abs(512.0-563.5)
ydist_712=abs(512.0-684.5)

##center for part 6 ROXS 42b
xdist_642=abs(512.0-626.0)
ydist_642=abs(512.0-484.0)

##center for part 6 ROXS 42b
xdist_612=abs(512.0-569.0)
ydist_612=abs(512.0-683.0)

##center for part 5 ROXS 42b
xdist_542=abs(512.0-626.0)
ydist_542=abs(512.0-477.0)

##center for part 5 ROXS 12
xdist_512=abs(512.0-570.0)
ydist_512=abs(512.0-683.0)

##center for part 4 ROXS 42
xdist_442=abs(512.0-628.5)
ydist_442=abs(512.0-513.0)

##center for part 4 ROXS 12
xdist_412=abs(512.0-483.0)
ydist_412=abs(512.0-690.0)

##redefine depending on which image i am looking at
xdist=xdist_742
ydist=ydist_742

dist_to_star= sqrt(xdist**2+ydist**2)
projected_sep = dist_to_star*arcsec_pixel

print (projected_sep)

# PA for ROXS 12
#theta_rel=arctan((xdist/ydist)) +pi/2.0

##PA for ROXS 42B
theta_rel=arctan(ydist/xdist)+pi

PA=theta_rel
print(degrees(PA))


figure(1)
plt.imshow(med_part7a,cmap='Greys_r')
plt.colorbar()
#
#figure(2)
#plt.imshow(rotated_median_12,cmap='Greys_r')
#plt.colorbar()

#figure(3)
#plt.imshow(ADI_med_12,cmap='Greys_r')
#plt.colorbar()
#
#figure(4)
#plt.imshow(ADI_sum_12,cmap='Greys_r')
#plt.colorbar()


show()
print(time.time()-starttime)
