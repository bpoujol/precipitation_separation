#ALGORITHM OF PRECIPITATION SEPARATION
#PART 2 : SEPARATION OF PRECIPITATION WITH THE THRESHOLDS
#This script can only be run once Part 1 has been completed.
#This step is fast but requires memory since many variables must be loaded at the same time.
#It is therefore recommended to use a subregion not too large if you work on the climatic time scales.

import numpy as np
from netCDF4 import Dataset

#Thresholds of the algorithm
#Threshold for vorticity
threshold_z = 0.5 #(m/s) (the threshold in s^-1 must be multiplied by the horizontal resolution of your data)
#Threshold for vertical velocity
threshold_w = 0.04 #(m/s)
#Threshold for topography-induced vertical velocity
threshold_topo = 0.2 #(m/s)

#Load the mask of the subregion where you want to separate the precipitation
mask = np.array(Dataset("mask_east.nc")["mask_array"])[:,:]
#Convert the mask into a boolean array
mask = mask_west.astype('bool')

#Load the precipitation data and mask it
precip_temp = np.array(Dataset("precip.nc")["RAINNC"])[:,:,:]
precip = np.array([x[mask_west] for x in precip_temp])
del precip_temp

#Load and mask the standard deviations produced by part 1
indexw = np.load("index_w.npy")
index_w = np.array([x[mask_west] for x in indexw])
del indexw
indexz = np.load("index_z.npy")
index_z = np.array([x[mask_west] for x in indexz])
del indexz
indextopo = np.load("index_topo.npy")
index_topo = np.array([x[mask_west] for x in indextopo])
del indextopo


#Separate the precipitation by creating boolean masks
#The logical expressions correspond to the figure 2 of the article

#CONVECTIVE PRECIPITATION
#Creation of the mask
mask_convective = np.logical_and(index_w>threshold_w,np.logical_or(index_topo<threshold_topo,np.logical_and(index_topo>threshold_topo,index_z>threshold_z)))
#This boolean expression corresponds to the paths that bring to convective precipitation on Fig2


#If your precipitation and your wind speed are not at the same time resolution, the next lines put them at the same resolution by duplicating the wind speed data.
#For example here precipitation is hourly and wind is k-hourly so all the masks must be duplicated three times
#Ensure that your input data files have corresponding time steps, i.e. the first wind speed time steps is in the middle of the k first time steps of the precipitation data.
k = 3
nt_dyn = np.shape(mask_convective)[0] #Number of timesteps in the wind speed data
#Resize the mask to the precipitation shape
indices = np.reshape(np.array([[i for j in range(k)]for i in range(nt_dyn)]),-1)
mask_convective = np.array([mask_convective[i] for i in indices])

#Separation of convective precipitation
#Make a copy of the precipitation data
convective = np.copy(precip)
#Put to zero all the non-convective precipitation rates
convective[np.logical_not(mask_convective)] = 0. 
#Delete a useless variable
del mask_convective
#Save the result
np.save("convective.npy",convective)
del convective


#OROGR-STRATI PRECIPITATION
#Same as convective precipitation.
mask_orographic = np.logical_and(index_w>threshold_w,np.logical_and(index_topo>threshold_topo,index_z<threshold_z))
mask_orographic = np.array([mask_orographic[i] for i in indices])
orographic = np.copy(precip)
orographic[np.logical_not(mask_orographic)] = 0.
del mask_orographic
np.save("orographic.npy",orographic)

#STRATIFORM PRECIPITATION
#Same as convective precipitation.
mask_stratiform = (index_w<threshold_w)
mask_stratiform = np.array([mask_stratiform[i] for i in indices])
stratiform = np.copy(precip)
stratiform[np.logical_not(mask_stratiform)] = 0.
del mask_stratiform
np.save("stratiform.npy",stratiform)
