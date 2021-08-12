# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 11:59:54 2021

@author: markl
Written by mark in 2021
Various functionalities of parcels and xarray to analyze parcels output
"""

from parcels import FieldSet, ParticleSet, JITParticle
from parcels import AdvectionRK4
import numpy as np
from datetime import timedelta as delta
import netCDF4
import xarray as xr
import matplotlib.pyplot as plt

##########SIMULATE PARTICLES AND EXPORT DATA TO NetCDF4###########
fieldset = FieldSet.from_parcels("parcels_examples/Peninsula_data/peninsula", allow_time_extrapolation=True)

npart = 10  # number of particles to be released
lon = 3e3 * np.ones(npart)
lat = np.linspace(3e3 , 45e3, npart, dtype=np.float32)
time = np.arange(0, npart) * delta(hours=2).total_seconds()  # release every particle two hours later

pset = ParticleSet(fieldset=fieldset, pclass=JITParticle, lon=lon, lat=lat, time=time)

output_file = pset.ParticleFile(name="Output.nc", outputdt=delta(hours=2))
pset.execute(AdvectionRK4, runtime=delta(hours=24), dt=delta(minutes=5),
             output_file=output_file)
#output_file.close()  # export the trajectory data to a netcdf file

######### READ IN NETCDF FILE THAT WAS CREATED ################
data_netcdf4 = netCDF4.Dataset('Output.nc')
print(data_netcdf4)

######### NAME VARIABLES FROM NETCDF FILE ################
trajectory_netcdf4 = data_netcdf4.variables['trajectory'][:]
time_netcdf4 = data_netcdf4.variables['time'][:]
lon_netcdf4 = data_netcdf4.variables['lon'][:]
lat_netcdf4 = data_netcdf4.variables['lat'][:]
print(trajectory_netcdf4)

####### OR USE XARRAY #################
data_xarray = xr.open_dataset('Output.nc')
print(data_xarray)
print(data_xarray['trajectory'])

###### TRAJECTORY DATA STRUCTURE #########
np.set_printoptions(linewidth=160)
ns_per_hour = np.timedelta64(1, 'h') # nanoseconds in an hour

print(data_xarray['time'].data/ns_per_hour) # time is stored in nanoseconds

##### ANALYSIS ##########
x = data_xarray['lon'].values
y = data_xarray['lat'].values
distance = np.cumsum(np.sqrt(np.square(np.diff(x))+np.square(np.diff(y))),axis=1)  # d = (dx^2 + dy^2)^(1/2)

real_time = data_xarray['time']/ns_per_hour # convert time to hours
time_since_release = (real_time.values.transpose() - real_time.values[:,0]) # substract the initial time from each timeseries

fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4),constrained_layout=True)

ax1.set_ylabel('Distance travelled [m]')
ax1.set_xlabel('observation',weight='bold')
d_plot = ax1.plot(distance.transpose())

ax2.set_ylabel('Distance travelled [m]')
ax2.set_xlabel('time since release [hours]',weight='bold')
d_plot_t = ax2.plot(time_since_release[1:],distance.transpose())
plt.show()

plt.figure()
ax= plt.axes()
ax.set_ylabel('Distance travelled [m]')
ax.set_xlabel('time [hours]',weight='bold')
d_plot_t = ax.plot(real_time.T[1:],distance.transpose())
