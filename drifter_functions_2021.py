# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 12:29:04 2021

@author: markl
A collection of functions to extract observational drifter data from
ERDDAP, plot the drifter on googlemaps and basemap. Originally made for BoF
Manuscript but generalizable.
"""

import pandas as pd
import sys
import numpy as np
from datetime import datetime as dt
import pytz
import urllib3
import pygmaps
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#To make inputs empty
#BoF_box = []
#time_range = []
#depth_range = []

#Your geographic, temporal, and depth ranges here
BoF_box = [-64.72, -67.35, 45.33, 44.23]  # maxlon, minlon,maxlat,minlat, for BoF range
time_range = [dt(2010, 11, 3, 0, 0, 0, 0, pytz.UTC),
              dt(2020, 12, 29, 0, 0, 0, 0, pytz.UTC)]  # start time and end time
depth_range = [-40, 40]  #min depth and max depth (meters)



def getobs_drift_by_info(gbox,input_time,depths):
    """
     Function written by Mark (modified from Huaxin's code)
     and used by plot_drift_gmap and plot_drift_basemap,(compatible with other functions)
     parameters: gbox (geographic box w/min and max lat/lon), input_time, and depths(min and max)
     return: ids, times, latitude,longitude, depths
     input_time can must contain two values: start_time & end_time
     example: input_time=[dt(2012,1,1,0,0,0,0,pytz.UTC),dt(2012,2,1,0,0,0,0,pytz.UTC)]
     Purpose: by user inputs (spatial, temporal, depth ranges), generate a URL 
     to find all drifter observations within specified ranges. Then return the
     relevant data from ERDDAP site from those obeservations
     """
    
    #Sequence of conditions to define ERDDAP URL based on all possible queries
    if gbox == [] and input_time == [] and depths == []:
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'
    elif gbox == [] and input_time == []:
        depth_min = depths[0]
        depth_max = depths[1]
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&time%3E='\
        +'&depth%3E='+str(depth_min) +'&depth%3C='+str(depth_max)
    elif gbox == [] and depths == []:
        time_min=input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
        time_max=input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')   
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&time%3E='+str(time_min)\
        +'&time%3C='+str(time_max)
    elif input_time == [] and depths == []:
        lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&latitude%3E='+str(lat_min)\
        +'&latitude%3C='+str(lat_max)+'&longitude%3E='+str(lon_min)\
        +'&longitude%3C='+str(lon_max)
    elif gbox == []:
        time_min=input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
        time_max=input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')   
        depth_min = depths[0]
        depth_max = depths[1]
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&time%3E='+str(time_min)\
        +'&time%3C='+str(time_max)+'&depth%3E='+str(depth_min)\
        +'&depth%3C='+str(depth_max)
    elif input_time == []:
        lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
        depth_min = depths[0]
        depth_max = depths[1]
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv0?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&latitude%3E='+str(lat_min)\
        +'&latitude%3C='+str(lat_max)+'&longitude%3E='+str(lon_min)\
        +'&longitude%3C='+str(lon_max)+'&depth%3E='+str(depth_min)\
        +'&depth%3C='+str(depth_max)
    elif depths == []:
        time_min=input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
        time_max=input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')   
        lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&time%3E='+str(time_min)\
        +'&time%3C='+str(time_max)+'&latitude%3E='+str(lat_min)\
        +'&latitude%3C='+str(lat_max)+'&longitude%3E='+str(lon_min)\
        +'&longitude%3C='+str(lon_max)
    else: 
        time_min=input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
        time_max=input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')   
        lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
        depth_min = depths[0]
        depth_max = depths[1]
        url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
        +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&time%3E='+str(time_min)\
        +'&time%3C='+str(time_max)+'&latitude%3E='+str(lat_min)\
        +'&latitude%3C='+str(lat_max)+'&longitude%3E='+str(lon_min)\
        +'&longitude%3C='+str(lon_max)+'&depth%3E='+str(depth_min)\
        +'&depth%3C='+str(depth_max)
    
            
    df = pd.read_csv(url,skiprows=[1]) #make drifter info table, USE THIS NORMALLY
    #df = pd.read_csv(r'C:\Users\markl\.spyder-py3\example_data.csv', skiprows=[1]) example for testing
    df.sort_values(by='time', ascending = True, inplace = True) #sort by date in ascending order
    
    #define outputs 
    ids = df.id.values
    times = df.time.values
    latitudes = df.latitude.values
    longitudes = df.longitude.values
    depths = df.depth.values
   
    #print(url) for testing
    return ids, times, latitudes, longitudes, depths



def getobs_drift_by_ids(identity):
    """
     Function written by Mark (modified from Huaxin's code) 
     and used by plot_drift_gmap and plot_drift_basemap (compatible with other functions)
     parameters: identity (drifter id#)
     return: ids, times, latitudes, longitudes, depths
     Purpose: based on user specified ID, output ALL drifter data from that deployment
     Not limited to geographic or temporal range
     """
    #Generate the URL
    url ='https://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csv?'\
    +'id%2Ctime%2Clatitude%2Clongitude%2Cdepth'+'&id%3E=%22'\
    +str(identity)+'%22&id%3C=%22'+str(identity)+'%22'


    df = pd.read_csv(url,skiprows=[1]) #drifter info table, USE THIS NORMALLY
    df.sort_values(by='time', ascending = True, inplace = True) #sort by date in ascending order
    
    #define outputs
    ids = df.id.values
    times = df.time.values
    latitudes = df.latitude.values
    longitudes = df.longitude.values
    depths = df.depth.values
    
    return ids, times, latitudes, longitudes,depths

def plot_drift_release_gmap():
    """
     Function written by Mark (modified from Huaxin's code)
     parameters: none
     output: saves a google map .html to your current file path, interactive online
     Purpose: Generate a google map of drifter RELEASES given that they pass through
     a certain geographic region, time, and depth range
     """
    ids, times, lats, lons, depths = getobs_drift_by_info(BoF_box,time_range,depth_range) #define inputs
    id = list(set(ids))
    #print(np.mean(lats)) for testing
    #print(np.mean(lons))
    print("There are " +str(len(id))+ " drifters that passed through the specified area and depth on their deployment.")
    
    mymap = pygmaps.pygmaps(np.mean(lats), np.mean(lons), 12) #initialize gmap
    
    for k in range(len(id)): #loop through the list of IDs generated by getobs_drift_by_info
    
        ids1, times1, lats1, lons1, depths1 = getobs_drift_by_ids(id[k]) #define new inputs by ID
        
        first_ping_time = np.min(times1) #earliest and latest pings
        last_ping_time = np.max(times1)
        first_ping_index = np.argmin(times1) #index of those pings
        last_ping_index = np.argmax(times1)
        first_ping_location = [lats1[first_ping_index],lons1[first_ping_index]] #position of those pings
        last_ping_location = [lats1[last_ping_index],lons1[last_ping_index]]
        
        #helpful info for the user 
        print("Drifter #" +str(id[k]) + " was released at " + str(first_ping_time)\
              + "(index " +str(first_ping_index)+ ") location:(" \
              + str(first_ping_location[0])+','+str(first_ping_location[1])\
              + ") and stopped transmitting at " + str(last_ping_time)\
              +"(index " + str(last_ping_index)+ ")  location:("\
              + str(last_ping_location[0])+','+str(last_ping_location[1]) + ")")
        
        for j in range(len(ids1)): #add point for the first ping of every deployment
            mymap.addpoint(lats1[first_ping_index], lons1[first_ping_index], 'black')
            
    #mymap.draw('./' + dt.datetime.now().strftime('%Y-%m-%d %H:%M') + '.html')
    mymap.draw('driftreleasemapALLTIME.html') #change this filename based on user
    
    

def plot_drift_tracks_gmap():
    """
    UPDATE THIS
     Function written by Mark (modified from Huaxin's code)
     parameters: none
     output: saves a google map .html to your current file path, interactive online
     Purpose: Generate a google map of drifter tracks given that they pass through
     a certain geographic region, time, and depth range
     """
    
    ids, times, lats, lons, depths = getobs_drift_by_info(BoF_box,time_range,depth_range) #define inputs
    id = list(set(ids))
    #print(np.mean(lats)) for testing
    #print(np.mean(lons))
    print("There are " +str(len(id))+ " drifters that passed through the specified area and depth on their deployment.")
    
    mymap = pygmaps.pygmaps(np.mean(lats), np.mean(lons), 12) #initialize gmap
    
    for k in range(len(id)): #loop through the list of IDs generated by getobs_drift_by_info
        #print(id[k]) for testing
        path = []
        ids1, times1, lats1, lons1, depths1 = getobs_drift_by_ids(id[k]) #define new inputs by ID
        for j in range(len(ids1)): #add point for each satellite ping
            #mymap.addpoint(lats1[0], lons1[0], 'black') this is for other function (release)
            mymap.addradpoint(lats1[j], lons1[j],295,'black')
            path.append((lats1[j], lons1[j]))
        mymap.addpath(path)
        #mymap.addradpoint(lats1[0], lons1[0], 295, "red") beginning
        #mymap.addradpoint(lats1[-1], lons1[-1], 295, "blue") end
        #
    #mymap.draw('./' + dt.datetime.now().strftime('%Y-%m-%d %H:%M') + '.html')
    mymap.draw('drifttrackmap.html') #change this filename based on user
    
    
    
    
def plot_drift_release_basemap():
    """
    UPDATE THIS -- to be fixed
     Function written by Mark 
     parameters: none
     output: basemap plot 
     Purpose: Generate a basemap of drifter RELEASES given that they pass through
     a certain geographic region, time, and depth range
     """
    ids, times, lats, lons, depths = getobs_drift_by_info(BoF_box,time_range,depth_range)
    id = list(set(ids))
    
    res = 'h'
    mymap = Basemap(projection='lcc', resolution=res, width = 2E5, height = 2E5, lat_0 =np.mean(lats), lon_0=np.mean(lons))
    mymap.fillcontinents(color = 'grey')
    mymap.drawcoastlines()
    mymap.drawmapboundary()
    
    for k in range(len(id)):
        #need to use the getobs_drift_by_id here to eliminate all other values.
        #FIX WHEN ERDDAP IS UP
        ids1, times1, lats1, lons1, depths1 = getobs_drift_by_ids(id[k])
        
        first_ping_time = np.min(times1) #earliest and latest pings
        last_ping_time = np.max(times1)
        first_ping_index = np.argmin(times1) #index of those pings
        last_ping_index = np.argmax(times1)
        first_ping_location = [lats1[first_ping_index],lons1[first_ping_index]] #position of those pings
        last_ping_location = [lats1[last_ping_index],lons1[last_ping_index]]
        
        #helpful info for the user 
        print("Drifter #" +str(id[k]) + " was released at " + str(first_ping_time)\
              + "(index " +str(first_ping_index)+ ") location:(" \
              + str(first_ping_location[0])+','+str(first_ping_location[1])\
              + ") and stopped transmitting at " + str(last_ping_time)\
              +"(index " + str(last_ping_index)+ ")  location:("\
              + str(last_ping_location[0])+','+str(last_ping_location[1]) + ")")
        
        for j in range(len(ids1)): #add point for the first ping of every deployment
            x,y = mymap(lons1[first_ping_index], lats1[first_ping_index])
            mymap.plot(x,y,'bo', markersize=3)
        
    #mymap.draw('./' + dt.datetime.now().strftime('%Y-%m-%d %H:%M') + '.html')
    plt.savefig('drifterreleasebasemap.png')
    plt.show()
    
    
    #MAKE A PLOT DRIFT TRACK BASEMAP PLOT_DRIFT_TRACK_BASEMAP
    

#def plot_drift_tracks_basemap():
    #
    
    
if __name__ == '__main__':
    # obs_id = 196440671 #drifter ID from ERDDAP
    # ids, times, obs_lat, obs_lon, depths = getobs_drift_by_ids(obs_id)
    # print(times)
    plot_drift_release_gmap()
    
    
    
    
    
    
    
    
    
    
    
    
    