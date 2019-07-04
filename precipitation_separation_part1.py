#ALGORITHM OF PRECIPITATION SEPARATION
#PART 1 : CALCULATION OF THE STANDARD DEVIATIONS
#This step requires most of the time

#Modules necessary for the algorithm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


#Load the data
#Your data must be 3D with the format precip[time,y,x] in mm/hr
#Load mid tropospheric vertical velocity
wa500 = np.array(Dataset("wa500.nc")["W_PL"])[:,0,:,:]

nt_dyn,nx,ny = np.shape(wa500) #Shape of the data


#PROCESS VERTICAL VELOCITY
#Initialization of the standard deviation matrix
index_w = np.zeros((nt_dyn,nx,ny))

#Half-length of the box in which the standard deviation is calculated. Must be expressed in number of grid points. Must be an integer.
r = 5

#Calculation of the standard deviation over each box
for i in range(r,nx-r) :
    for j in range(r,ny-r):
        w_temp = wa500[:,i-r:i+r,j-r:j+r]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))
#Special cases at the borders of the domain
for i in range(r) :
    for j in range(r,ny-r):
        w_temp = wa500[:,0:i+r,j-r:j+r]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))
for i in range(nx-r,nx) :
    for j in range(r,ny-r):
        w_temp = wa500[:,i-r:nx,j-r:j+r]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r,nx-r) :
    for j in range(r):
        w_temp = wa500[:,i-r:i+r,0:j+r]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r,nx-r) :
    for j in range(ny-r,ny):
        w_temp = wa500[:,i-r:i+r,j-r:ny]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r) :
    for j in range(r):
        w_temp = wa500[:,0:i+r,0:j+r]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r) :
    for j in range(ny-r,ny):
        w_temp = wa500[:,0:i+r,j-r:ny]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(r,ny-r):
        w_temp = wa500[:,i-r:nx,0:j+r]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(ny-r,ny):
        w_temp = wa500[:,i-r:nx,j-r:ny]
        index_w[:,i,j] = np.std(w_temp,axis=(1,2))

#Delete useless varibles and save the standard deviation.
del wa500
np.save("index_w.npy",index_w)
del index_w

#PROCESS TOPOGRAPHY

#Load the topography. Must be a 2D array on the same grid than the wind speed data.
topo = np.array(Dataset("topography.nc")["Band1"])
#Load the 700 hPa wind speed along the x coordinate of the grid
ua700 = np.array(Dataset("u700.nc")["U_PL"])[:,0,:,:]
#Contribution of the zonal wind speed to the uplift vertical velocity
wtopo = ua700*np.array(nt_dyn*[np.gradient(topo,axis=1)])
#Delete the wind speed along x, now useless
del ua700
#Load the 700hPa wind speed along the y coordinate of the grid
va700 = np.array(Dataset("v700.nc")["V_PL"])[:,0,:,:]
#Add the contribution of meridional wind speed to the uplift velocity
wtopo += va700*np.array(nt_dyn*[np.gradient(topo,axis=0)])
del va700
del topo

#Divide this value by the horizontal resolution of your data (in meters) to convert it into m/s
wtopo = wtopo/3000.

#Initialization of the standard deviations
index_topo = np.zeros((nt_dyn,nx,ny))

#Calculation of the standard deviations
for i in range(r,nx-r) :
    for j in range(r,ny-r):
        w_temp = wtopo[:,i-r:i+r,j-r:j+r]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r) :
    for j in range(r,ny-r):
        w_temp = wtopo[:,0:i+r,j-r:j+r]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(r,ny-r):
        w_temp = wtopo[:,i-r:nx,j-r:j+r]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r,nx-r) :
    for j in range(r):
        w_temp = wtopo[:,i-r:i+r,0:j+r]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r,nx-r) :
    for j in range(ny-r,ny):
        w_temp = wtopo[:,i-r:i+r,j-r:ny]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r) :
    for j in range(r):
        w_temp = wtopo[:,0:i+r,0:j+r]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(r) :
    for j in range(ny-r,ny):
        w_temp = wtopo[:,0:i+r,j-r:ny]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(r,ny-r):
        w_temp = wtopo[:,i-r:nx,0:j+r]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(ny-r,ny):
        w_temp = wtopo[:,i-r:nx,j-r:ny]
        index_topo[:,i,j] = np.std(w_temp,axis=(1,2))

#Delete useless variables and save the results
del wtopo
np.save("index_topo.npy",index_topo)
del index_topo

#PROCESS VORTICITY
#Load 500 hPa wind speed along the x direction
ua500 = np.array(Dataset("u500.nc")["U_PL"])[:,0,:,:]
#Contribution of the zonal wind to relative vorticity
zeta =  -np.gradient(ua500,axis=1)
#Delete the zonal wind
del ua500

#Load 500 hPa wind along the y direction
va500 = np.array(Dataset("v500.nc")["V_PL"])[:,0,:,:]
#Add contribution of the meridional wind to relative vorticity
zeta += np.gradient(va500,axis=2)
#Delete the meridional wind
del va500

#Initialize and calculate the standard deviations
index_z = np.zeros((nt_dyn,nx,ny))
for i in range(r,nx-r) :
    for j in range(r,ny-r):
        z_temp = zeta[:,i-r:i+r,j-r:j+r]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))
for i in range(r) :
    for j in range(r,ny-r):
        z_temp = zeta[:,0:i+r,j-r:j+r]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))
for i in range(nx-r,nx) :
    for j in range(r,ny-r):
        z_temp = zeta[:,i-r:nx,j-r:j+r]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

for i in range(r,nx-r) :
    for j in range(r):
        z_temp = zeta[:,i-r:i+r,0:j+r]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

for i in range(r,nx-r) :
    for j in range(ny-r,ny):
        z_temp = zeta[:,i-r:i+r,j-r:ny]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

for i in range(r) :
    for j in range(r):
        z_temp = zeta[:,0:i+r,0:j+r]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

for i in range(r) :
    for j in range(ny-r,ny):
        z_temp = zeta[:,0:i+r,j-r:ny]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(r,ny-r):
        z_temp = zeta[:,i-r:nx,0:j+r]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

for i in range(nx-r,nx) :
    for j in range(ny-r,ny):
        z_temp = zeta[:,i-r:nx,j-r:ny]
        index_z[:,i,j] = np.std(z_temp,axis=(1,2))

del zeta

#Save the standard deviation
np.save("index_z.npy",index_z)
