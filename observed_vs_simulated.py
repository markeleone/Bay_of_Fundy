# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 12:51:09 2021

@author: markl
Written by mark in 2021
All functions to compare simulated and observed drifters and make graphs of the 
error and discrepancies
"""
import numpy as np
import matplotlib.tri as Tri 
import matplotlib.path as mpltpath
from datetime import datetime as dt
from datetime import timedelta as delta
import pandas as pd
import netCDF4 as nc
import os
import pathlib as path
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj") #needed to run basemap
from mpl_toolkits.basemap import Basemap
from scipy import stats
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile, AdvectionAnalytical, ScipyParticle
from drifter_functions_2021 import getobs_drift_by_ids
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import haversine as hs

##### HARDCODES FOR USER #########
obs_id = 196440671 #drifter ID from ERDDAP
ids, times, obs_lat, obs_lon, depths = getobs_drift_by_ids(obs_id)
print(times)





######MAIN CODE FOR DOPPIO SIMULATIONS#################



filenames = {'U': 'DOPPIO_2019_06_16_2019_07_07.nc',
              'V': 'DOPPIO_2019_06_16_2019_07_07.nc'}
variables = {'U': 'u_eastward',
              'V': 'v_northward'}
dimensions = {'lat':'lat_rho',
              'lon':'lon_rho',
              'time':'time1'}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
"""GET RUNGE-KUTTA 4 PSET"""



data_xarray = xr.open_dataset('BoFParticlesRK4_2.nc')
print(data_xarray)
print(data_xarray['trajectory'])

#np.set_printoptions(linewidth=160)
#ns_per_hour = np.timedelta64(1, 'h') # nanoseconds in an hour

#print('time datatype ' + str(type(data_xarray['time'])))
#print('operand data type ' + str(type(ns_per_hour)))

#print(data_xarray['time'].data/ns_per_hour) # time is stored in nanoseconds

x = data_xarray['lon'].values
y = data_xarray['lat'].values
distance = np.cumsum(np.sqrt(np.square(np.diff(x))+np.square(np.diff(y))),axis=1)  # d = (dx^2 + dy^2)^(1/2)


#real_time = data_xarray['time'].dt.day #gives day of the month, could be useful
real_time = data_xarray['time']#.astype('timedelta64[h]')


#real_time = data_xarray['time']/ns_per_hour # convert time to second
time_since_release = (real_time.values.transpose() - real_time.values[:,0]) # substract the initial time from each timeseries USE THIS WITH MULTIPLE RELEASE TIMES


# fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4),constrained_layout=True)

# """To plot observation # vs distance traveled"""
# ax1.set_ylabel('Distance travelled [m]')
# ax1.set_xlabel('observation',weight='bold')
# d_plot = ax1.plot(distance.transpose())

# """To plot time since release vs distance traveled on same figure"""
# ax2.set_ylabel('Distance travelled [m]')
# ax2.set_xlabel('time since release [day]',weight='bold')
# d_plot_t = ax2.plot(time_since_release[1:],distance.transpose()) #USE THIS WITH MULTIPLE RELEASE TIMES


# """To plot time (overall) vs distance traveled on new figure"""
# plt.figure()
# ax= plt.axes()
# ax.set_ylabel('Distance travelled [m]')
# ax.set_xlabel('time [hours]',weight='bold')
# d_plot_t = ax.plot(real_time.T[1:],distance.transpose())



# """To plot actual points along the lat-long axis"""
# fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(16,3.5),constrained_layout=True)

# ###-Points-###
# ax1.set_title('Points')
# ax1.scatter(data_xarray['lon'].T,data_xarray['lat'].T)
# ###-Lines-###
# ax2.set_title('Lines')
# ax2.plot(data_xarray['lon'].T,data_xarray['lat'].T)
# ###-Points + Lines-###
# ax3.set_title('Points + Lines')
# ax3.plot(data_xarray['lon'].T,data_xarray['lat'].T,marker='o')
# ###-Not Transposed-###
# ax4.set_title('Not transposed')
# ax4.plot(data_xarray['lon'],data_xarray['lat'],marker='o')


"""To plot observed vs simulated drifter"""
#read in lat and lon
lat = data_xarray['lat'].values
lon = data_xarray['lon'].values

lat_min_sim = min(min(lat))
lat_max_sim = max(max(lat))
lon_min_sim = min(min(lon))
lon_max_sim = max(max(lon))
lat_min_obs = min(obs_lat)
lat_max_obs = max(obs_lat)
lon_min_obs = min(obs_lon)
lon_max_obs = max(obs_lon)

# determine range to print based on min, max lat and lon of the data
margin = 0.2 # buffer to add to the range
lat_min = min(lat_min_sim,lat_min_obs) - margin
lat_max = max(lat_max_sim, lat_max_obs) + margin
lon_min = min(lon_min_sim, lon_min_obs) - margin
lon_max = max(lon_max_sim, lon_max_obs) + margin


# create map using BASEMAP
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            lat_0=(lat_max - lat_min)/2,
            lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'h',
            area_thresh=10000.,
            )
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.drawmapboundary(fill_color='#46bcec')
m.fillcontinents(color = 'white',lake_color='#46bcec')
# convert lat and lon to map projection coordinates
lons, lats = m(lon, lat)
obs_lons, obs_lats = m(obs_lon, obs_lat)
# plot points as red dots
m.plot(lons.T, lats.T, color='r', zorder=6, label='Simulated (DOPPIO)') #add marker 'o' to get dots, .T connects line #FOR SIMULATED
m.plot(obs_lons.T, obs_lats.T, color='b', zorder=5, label='Observed') #OBSERVED
plt.legend(loc='lower right', shadow=True)
plt.title('June 16 to July 7 2019')

m.drawparallels(np.arange(-90,90,0.2),labels=[1,1,0,1], fontsize=8)
m.drawmeridians(np.arange(-180,180,0.2),labels=[1,1,0,1], rotation=45, fontsize=8)
# draw parallels and meridians.
# label parallels on right and top
# meridians on bottom and left

plt.show()

print("Simulated observations: " + str(lon.size))
print("Observed observations: " + str(len(obs_lon)))

#plotTrajectoriesFile('BoFParticlesRK4.nc', mode='movie2d_notebook') #recieved error message
#plt = plotTrajectoriesFile('BoFParticlesRK4.nc')


                      # tracerfile='his_History_Best.nc',
                      # tracerlon='lon_rho',
                      # tracerlat='lat_rho',
                      # tracerfield='ubar_eastward')"""



# particle_trajectories = 'BoFParticles.nc'
# ds = nc.Dataset(particle_trajectories)

# print(ds)

def separation_distance(obs_lat,obs_lon,sim_lat,sim_lon):
    """Function written by Mark in 2021 to create a graph of separation distance
    (not normalized by trajectory length) over time
    input variables are numpy arrays of any length representing 
    observed and simulated latitude and longitude"""
    obs_lat = obs_lat.tolist(); obs_lon = obs_lon.tolist()
    sim_lat = sim_lat.flatten().tolist(); sim_lon = sim_lon.flatten().tolist()
    #truncate observations to 72 hr
    list_size = 72
    obs_lat = obs_lat[0:list_size]; obs_lon = obs_lon[0:list_size]
    sim_lat = sim_lat[0:list_size]; sim_lon = sim_lon[0:list_size]
    print(len(sim_lat))
    
    hs_distance = []
    for i in range(len(obs_lat)):
        obs_loc = (obs_lat[i], obs_lon[i])
        sim_loc = (sim_lat[i], sim_lon[i])
        hs_distance.append(hs.haversine(obs_loc,sim_loc))#default returns distances in kilometers
        
    print(hs_distance) 
    x=np.linspace(0,72,72)
    plt.plot(x,hs_distance,'o')
    plt.title('Separation distance over time')
    plt.xlabel('time (hours)')
    plt.ylabel('separation distance (km)')
    
    
    
if __name__ == '__main__':
    separation_distance(obs_lat,obs_lon,lat,lon)

