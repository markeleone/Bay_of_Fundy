# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 13:07:04 2021

@author: markl
A collection of functions for visualizing FVCOM velocity fields, with tidal inputs and detided.
Originally made for the GOM4 grid for the BoF manuscript by Mark. Code also contributed
by JiM and Vitalii as credited in function descriptions.
"""

from pylab import *
import matplotlib.tri as Tri
from datetime import datetime as dt
import pandas as pd
import pytz
import glob
import numpy as np
import netCDF4
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj") #needed to run basemap
from mpl_toolkits.basemap import Basemap


def inconvexpolygon(xp,yp,xv,yv):
    """
	check if point is inside a convex polygon

    i=inconvexpolygon(xp,yp,xv,yv)
    
    xp,yp - arrays of points to be tested
    xv,yv - vertices of the convex polygon

    i - boolean, True if xp,yp inside the polygon, False otherwise
    From JiM's uvmodel.py
    """    
    N=len(xv)   
    j=np.arange(N)
    ja=(j+1)%N # next vertex in the sequence    
    NP=len(xp)
    i=np.zeros(NP,dtype=bool)
    for k in range(NP):
        # area of triangle p,j,j+1
        Aj=np.cross(np.array([xv[j]-xp[k],yv[j]-yp[k]]).T,np.array([xv[ja]-xp[k],yv[ja]-yp[k]]).T) 
    	# if a point is inside the convect polygon all these Areas should be positive 
    	# (assuming the area of polygon is positive, counterclockwise contour)
        Aj /= Aj.sum()
    	# Now there should be no negative Aj
    	# unless the point is outside the triangular mesh
        i[k]=(Aj>0.).all()
        
    return i


def nearlonlat(lon,lat,lonp,latp): # needed for the next function get_FVCOM_bottom_temp
    """
    i=nearlonlat(lon,lat,lonp,latp) change
    find the closest node in the array (lon,lat) to a point (lonp,latp)
    input:
        lon,lat - np.arrays of the grid nodes, spherical coordinates, degrees
        lonp,latp - point on a sphere
        output:
            i - index of the closest node
            For coordinates on a plane use function nearxy           
            from Vitalii Sheremet, FATE Project via JiM's uvmodel.py
    """
    cp=np.cos(latp*np.pi/180.)
    # approximation for small distance
    dx=(lon-lonp)*cp
    dy=lat-latp
    dist2=dx*dx+dy*dy
    i=np.argmin(dist2)
    return i


def get_FVCOM_data(lati,loni,dtime,layer): # gets modeled temp using GOM3 forecast
        '''
        Taken primarily from Rich's blog at: http://rsignell-usgs.github.io/blog/blog/2014/01/08/fvcom/ on July 30, 2018
        where lati and loni are the position of interest, dtime is the datetime, and layer is "-1" for bottom
        '''
        #urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
        #urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
        urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM4_FORECAST.nc' #use GOM4 (Mark)
        nc = netCDF4.Dataset(urlfvcom).variables
        #first find the index of the grid 
        lat = nc['lat'][:]
        lon = nc['lon'][:]
        #second find the index of time
        time_var = nc['time']
        itime = netCDF4.date2index(dtime,time_var,select='nearest')# where startime in datetime
        inode = nearlonlat(lon,lat,loni,lati)
        temp = nc['temp'][itime,layer,inode]
        u = nc['u'][itime,layer,inode] #eastward water velocity (m/s)
        v = nc['v'][itime,layer,inode] #northward water velocity (m/s)
        #ww is upward water velocity
        #h is sea floor depth below geoid (bathymetry) units: m
        return temp, u, v
    
    
def display_FVCOM_data():
    "Written by Mark for understanding thermocline at a given position and a time"
    #set parameters
    latitude = 44.32
    longitude = -67.83
    datetime = dt(2021, 6, 22, 0, 0, 0, 0)
    bottom = -1
    surface = 0
    #Just to check for stratification in cape cod bay to make sure bottom layer actually = -1. It does
    bottom_temp, bottom_u, bottom_v = get_FVCOM_data(latitude, longitude, datetime, bottom)
    surface_temp, surface_u, surface_v = get_FVCOM_data(latitude, longitude, datetime, surface)
    print('Surface temp: ' + str(surface_temp) + ' ------ Bottom temp: ' + str(bottom_temp))
    print('Surface u: ' + str(surface_u) + ' ------ Bottom u: ' + str(bottom_u))
    print('Surface v: ' + str(surface_v) + ' ------ Bottom v: ' + str(bottom_v))
    
    
def sh_bindata(x, y, z, xbins, ybins):
    """ From uvmodel.py
    for binning velocity field vector data (for both u and v components of velocity)
    given xyz position and the x and y bins, return the location (xb,yb)
    and the new component of velocity data (zb_mean)
    """
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
    return xb,yb,zb_mean


def sh_rmtide(f,dt=1.,ends=0.):
    """
    removes solar and lunar tidal signal by sequentially averaging 
    over their periods: of 24h and 24h50.47m. This is equivalent to
    applying convolution with a trapezoidal shaped window (almost triagular).
    
    f - input tidal sequence
    dt - uniform time step of input series, in hours, dt=1. hr default
    ends - fill value for the ends: 0. or nan
    
    fm = sh_rmtide(f,dt=1.,ends=np.nan)
    from uvmodel.py
    """
    TS=24. # Solar hour angle period 24h
    TL=24.+50.47/60. # Lunar hour angle period
    N=len(f)
    fm=np.zeros(N)*ends

    # remove solar period    
    T=TS/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
    # all weights =1 except at ends      
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
        
    # remove lunar period
    f=fm*1.0  # deep copy!  
    T=TL/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
    # all weights =1 except at ends      
    #    print w
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
    
    return fm


def get_FVCOM_GOM4_url(month,year):
    """Give GOM4 url for a specific month by Mark"""
    url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM/' + year + '/gom4_' + year + month +'.nc'
    return url
    

def plot_FVCOM_velocity_field():
    """Most of code inspired by Rich Signell's blog at http://rsignell-usgs.github.io/blog/blog/2014/01/08/fvcom/
    Modified by Mark in summer 2021
    This function plots FVCOM velocity vector field in a given area 
    Optional to plot depth as another layer shaded (colorbar)"""
    # (rich's) url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
    #url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM4_FORECAST.nc' #forecast
    url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM/2019/gom4_201905.nc' #2019 data
    nc = netCDF4.Dataset(url).variables
    
    # Desired time for snapshot
    # ....right now (or some number of hours from now) ...
    start = dt(2019, 5, 15, 0, 0) #+ dt.delta(hours=18)
    # ... or specific time (UTC)
    #start = dt.datetime(2013,3,2,15,0,0)
    # Get desired time step  
    time_var = nc['time']
    itime = netCDF4.date2index(start,time_var,select='nearest')
    
    # Get lon,lat coordinates for nodes (depth)
    lat = nc['lat'][:]
    lon = nc['lon'][:]
    # Get lon,lat coordinates for cell centers (depth)
    latc = nc['latc'][:]
    lonc = nc['lonc'][:]
    # Get Connectivity array
    nv = nc['nv'][:].T - 1 
    # Get depth
    h = nc['h'][:]  # depth
    
    dtime = netCDF4.num2date(time_var[itime],time_var.units)
    daystr = dtime.strftime('%Y-%b-%d %H:%M')

    tri = Tri.Triangulation(lon,lat, triangles=nv)
    
    # get current at layer [0 = surface, -1 = bottom]
    ilayer = 0
    u = nc['u'][itime, ilayer, :]
    v = nc['v'][itime, ilayer, :]
    
    #bof
    levels=arange(-250,2,1)   # depth contours to plot
    ax = [-67.35, -64.72, 44.23, 45.33]
    
    #ax= [-70.97, -70.82, 42.25, 42.35] #boston harbor rich's
    maxvel = 0.5
    subsample = 6
    
    # find velocity points in bounding box
    ind = argwhere((lonc >= ax[0]) & (lonc <= ax[1]) & (latc >= ax[2]) & (latc <= ax[3]))
    
    np.random.shuffle(ind)
    Nvec = int(len(ind) / subsample)
    idv = ind[:Nvec]

    # tricontourf plot of water depth with vectors on top
    figure(figsize=(18,10))
    subplot(111,aspect=(1.0/cos(mean(lat)*pi/180.0)))
    tricontourf(tri, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
    axis(ax)
    gca().patch.set_facecolor('0.5')
    cbar=colorbar()
    cbar.set_label('Water Depth (m)', rotation=-90)
    Q = quiver(lonc[idv],latc[idv],u[idv],v[idv],scale=20)
    maxstr='%3.1f m/s' % maxvel
    qk = quiverkey(Q,0.92,0.08,maxvel,maxstr,labelpos='W')
    title('NECOFS Velocity, Layer %d, %s UTC' % (ilayer, daystr))
    plt.savefig("NECOFS_Velocity_depth-1_05_15_2019.png") #edit when changing time or space


def plot_FVCOM_binned_velocity_field():
    """Adapted from uvmodel.py updated by Mark in summer 2021
    plots FVCOM velocity field binned in time and space by calling function sh_bindata
    The output is a quiver plot that is tidally variable"""
    #urlname= 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM4_FORECAST.nc' #forecast/nowcast
    #urlname ='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM/2019/gom4_201905.nc' #2019
    month = '05'
    year = '2016'
    urlname = get_FVCOM_GOM4_url(month, year)
    #depth=-1 #surface
    depth=-1
    START_TIME= '2016-05-01 00:00:00' #CHANGE THIS
    END_TIME = '2016-05-15 00:00:00'
    START_TIME=dt.strptime(START_TIME, "%Y-%m-%d %H:%M:%S") 
    END_TIME=dt.strptime(END_TIME, "%Y-%m-%d %H:%M:%S") 
    # START_TIME = dt(2019, 5, 2, 0, 0)
    # END_TIME = dt(2019, 5, 17, 0, 0)
    timeperiod = (END_TIME)-(START_TIME)
    startrecord=int(((timeperiod.total_seconds())/60/60))
    

    print("timeperiod= "+str(timeperiod))
    print("startrecord= "+str(startrecord))
    lon_bin_size= 0.05
    lat_bin_size= 0.05
    gbox=[-67.5, -65, 44, 45.3] #Bay of Fundy
    
    #plot_FVCOM_velocity_field() this works
    print(urlname)
    nc = netCDF4.Dataset(urlname)
    print("hold on readin nc")
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    latc = nc.variables['latc'][:]
    lonc = nc.variables['lonc'][:]
    h = nc.variables['h'][:]
    siglay=nc.variables['siglay']
    u= nc.variables['u']
    v= nc.variables['v']
    nv = nc.variables['nv'][:].T - 1
    print(u)

        
    print('we have the model data now we want to get depth of interest')
    utotal=[]
    vtotal=[]
    if depth==-1: # case of surface flow
      utotal=u[startrecord,0,:]
      vtotal=v[startrecord,0,:]
      print(utotal)
      print(vtotal)
    else:
      for i in range(len(lon)):
        depthtotal=siglay[:,i]*h[i]
        layer=np.argmin(abs(depthtotal+depth))
        utotal.append(u[startrecord,layer,i])
        vtotal.append(v[startrecord,layer,i])
      utotal=np.array(utotal)
      vtotal=np.array(vtotal)
    
    print('now lets bin the data')
    #xi = np.arange(min(lon)-0.1,max(lon)+0.1,lon_bin_size)
    #yi = np.arange(min(lat)-0.1,max(lat)+0.1,lat_bin_size)
    xi = np.arange(gbox[0]-0.1,gbox[1]+0.1,lon_bin_size)
    yi = np.arange(gbox[2]-0.1,gbox[3]+0.1,lat_bin_size)
    xb,yb,ub_mean = sh_bindata(lon[::-1], lat[::-1], utotal, xi, yi)
    xb,yb,vb_mean = sh_bindata(lon[::-1], lat[::-1], vtotal, xi, yi)
    xxb,yyb = np.meshgrid(xb, yb)
    #latsize=[min(lat)-0.6,max(lat)+0.6]
    #lonsize=[min(lon)-0.6,max(lon)+0.6]
    latsize=[gbox[2]-0.1,gbox[3]+0.1]
    lonsize=[gbox[0]-0.1,gbox[1]+0.1]
    print('and plot')
    plt.figure()
    m = Basemap(projection='cyl',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
                llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='h')#,fix_aspect=False)
    #m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,5),labels=[1,0,0,0])
    #m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,10),labels=[0,0,0,1])
    m.drawparallels(np.arange(min(latsize),max(latsize),10),labels=[1,0,0,0])
    m.drawmeridians(np.arange(min(lonsize),max(lonsize),10),labels=[0,0,0,1])
    m.drawcoastlines()
    m.fillcontinents(color='grey')
    m.drawmapboundary()
    ub = np.ma.array(ub_mean, mask=np.isnan(ub_mean))
    vb = np.ma.array(vb_mean, mask=np.isnan(vb_mean))
    Q=m.quiver(xxb,yyb,ub,vb,scale=20)
    
    plt.quiverkey(Q,0.8,0.05,1, '1m/s', labelpos='W')
    plt.title('FVCOM Velocity Field Depth:'+str(depth)+' Time:'+str(START_TIME)) 
    plt.savefig('FVCOM_velocity_field_with_tides_05022016-05152016_15m.png')
    plt.show()
    
    
    
def plot_FVCOM_binned_detided_velocity_field():
    """Adapted from uvmodel.py updated by Mark in summer 2021
    plots FVCOM velocity field binned in time and space by calling function sh_bindata
    removes tidal influence by calling function sh_rmtide
    The output is a quiver plot that is not tidally variable"""
    #urlname= 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM4_FORECAST.nc' #forecast
    #urlname ='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM/2019/gom4_201905.nc' #2019
    month = '05'
    year = '2016'
    urlname = get_FVCOM_GOM4_url(month, year)
    #depth=-1 #surface
    depth=-1
    START_TIME= '2016-05-01 00:00:00' #CHANGE THIS
    END_TIME = '2016-05-15 00:00:00'
    START_TIME=dt.strptime(START_TIME, "%Y-%m-%d %H:%M:%S") 
    END_TIME=dt.strptime(END_TIME, "%Y-%m-%d %H:%M:%S") 
    # START_TIME = dt(2019, 5, 2, 0, 0)
    # END_TIME = dt(2019, 5, 17, 0, 0)
    timeperiod = (END_TIME)-(START_TIME)
    startrecord=int(((timeperiod.total_seconds())/60/60))
    
    print("timeperiod= "+str(timeperiod))
    print("startrecord= "+str(startrecord))
    lon_bin_size= 0.05
    lat_bin_size= 0.05
    gbox=[-67.5, -65, 44, 45.3] #Bay of Fundy
    detide = 'yes'
        
    num_days=1 #number of days to make picture
    
    print(urlname)
    nc = netCDF4.Dataset(urlname)
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    latc = nc.variables['latc'][:]
    lonc = nc.variables['lonc'][:]
    h = nc.variables['h'][:]
    siglay=nc.variables['siglay']
    u= nc.variables['u']
    v= nc.variables['v']
    nv = nc.variables['nv'][:].T - 1
    print('we have the model data and now we are finding indices within the gbox')
    gbox_poly=array([[gbox[0],gbox[2]],[gbox[1],gbox[2]],[gbox[1],gbox[3]],[gbox[0],gbox[3]]])
    i=inconvexpolygon(lonc,latc,gbox_poly[:,0],gbox_poly[:,1])
    i=np.argwhere(i).flatten()
    llond0=lonc[i] # reduces lat/lon to inside the polygon
    llatd0=latc[i]
    num_model_pts_in_gbox=len(i)
    print('number of model points inside gbox = '+str(num_model_pts_in_gbox))
    #THIS ONLY CODED FOR THE DETIDED CASE!!!!!!!!!!!!!!!!!!!!!!!
    for n in range(num_days): # though through one day at a time (increments by 24 hours at end of loop)
        utotal=[0]*num_model_pts_in_gbox
        vtotal=[0]*num_model_pts_in_gbox
        utotal=np.zeros(num_model_pts_in_gbox,dtype=float)
        vtotal=np.zeros(num_model_pts_in_gbox,dtype=float)
        
        if depth==-1: # case of surface flow
            if detide=='no':
                utotal=u[startrecord,0,:] 
                vtotal=v[startrecord,0,:] 
                
            else: 
                print('detiding')
    	    	#for k in range(len(u[startrecord,0,np.argwhere(i).flatten()])): # for each velocity inside the polygon, detide
                for k in range(num_model_pts_in_gbox): # for each velocity inside the polygon, detide
                    utotal[k]=sh_rmtide(u[startrecord-24:startrecord+24,0,i[k]])[24]# uses 48 hours to apply a filter around both solar and lunar periods
                    vtotal[k]=sh_rmtide(v[startrecord-24:startrecord+24,0,i[k]])[24] 
        else:
            for i in range(len(lon)): 
                depthtotal=siglay[:,i]*h[i]
                layer=np.argmin(abs(depthtotal+depth)) 
                utotal=np.append(utotal,u[startrecord,layer,i]) 
                vtotal=np.append(vtotal,v[startrecord,layer,i])
        utotal=np.array(utotal) 
        vtotal=np.array(vtotal)
        print(utotal)
        print(vtotal)
        
        print('now lets bin the data')
        xi = np.arange(gbox[0]-0.1,gbox[1]+0.1,lon_bin_size)
        yi = np.arange(gbox[2]-0.1,gbox[3]+0.1,lat_bin_size)
        xb,yb,ub_mean = sh_bindata(llond0[::-1], llatd0[::-1], utotal, xi, yi)
        xb,yb,vb_mean = sh_bindata(llond0[::-1], llatd0[::-1], vtotal, xi, yi) 
        xxb,yyb = np.meshgrid(xb, yb)
        latsize=[gbox[2]-0.1,gbox[3]+0.1] 
        lonsize=[gbox[0]-0.1,gbox[1]+0.1] 
        
        print('and plot') 
        plt.figure() 
        m = Basemap(projection='cyl',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
                llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='i')#,fix_aspect=False) 
        m.drawparallels(np.arange(int(min(latsize)),int(ceil(max(latsize))),round(diff(latsize)[0]/4,1)),labels=[1,0,0,0]) 
        m.drawmeridians(np.arange(int(min(lonsize)),int(ceil(max(lonsize))),round(diff(lonsize)[0]/4,1)),labels=[0,0,0,1])
    	#m.drawparallels(np.arange(min(latsize),max(latsize),5),labels=[1,0,0,0])
    	#m.drawmeridians(np.arange(min(lonsize),max(lonsize),10),labels=[0,0,0,1]) 
        m.drawcoastlines() 
        m.fillcontinents(color='grey') 
        m.drawmapboundary() 
        ub = np.ma.array(ub_mean, mask=np.isnan(ub_mean)) 
        vb = np.ma.array(vb_mean, mask=np.isnan(vb_mean)) 
        Q=m.quiver(xxb,yyb,ub,vb,scale=5) 
        
        plt.quiverkey(Q,0.8,0.05,1, '1 m/s', labelpos='W') 
        plt.title('FVCOM velocity field detided')
        plt.savefig('FVCOM_velocity_field_detided_05012016-05152016_30m.png')
        plt.show() 
        print('Saved figure as <grid>_vectors<TIME>') 
        
        #plt.close('all') 
        startrecord=startrecord+24# jumps one day
    



if __name__ == '__main__':
    plot_FVCOM_binned_velocity_field()
    