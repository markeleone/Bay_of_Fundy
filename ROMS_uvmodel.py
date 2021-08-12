# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:03:01 2021

@author: markl

Inspired by PICO github, uvmodel.py and Rich Signell's blog
Still working on binned and detided data processing
"""

import numpy as np
from pylab import *
import matplotlib.tri as Tri 
import matplotlib.path as mpltpath
from datetime import datetime as dt
from datetime import timedelta
from dateutil import parser
import pandas as pd
import pytz
import glob
import os
import pathlib as path
from mpl_toolkits.basemap import Basemap
from scipy import stats
import urllib3
import netCDF4


#url = 'http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd' #COAWST
#url = 'https://tds.marine.rutgers.edu/thredds/dodsC/floats/doppio_flt_20200731.nc' #DOPPIO
#url = 'http://testbedapps-dev.sura.org/thredds/dodsC/alldata/Shelf_Hypoxia/tamu/roms/tamu_roms.nc'
#url = 'http://tds.ve.ismar.cnr.it:8080/thredds/dodsC/ismar/model/field2/run1/Field_2km_1108_out30min_his_0724.nc'
#url='http://tds.ve.ismar.cnr.it:8080/thredds/dodsC/field2_test/run1/his'

###HARDCODES FOR getroms_field_byinfo()
BoF_box = [-64.7255, -67.3578, 45.3302, 44.2376]  # maxlon, minlon,maxlat,minlat, for BoF range
time_range =  ['2019', '06', '16', '01', '00', '00',
               '2019', '06', '19', '12', '00', '00']  # start time and end time
s_rho = -0.9875  #-1 is case of surface flow


def getroms_field_byinfo(gbox,input_time,s_rho):
    """to fetch correct doppio url and velcoity field for an actual time range 
    updated version of det_doppio_url written by Mark in 2021"""
    lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
    min_year = input_time[0]; max_year = input_time[6]
    min_month = input_time[1]; max_month = input_time[7]
    min_day = input_time[2]; max_day = input_time[8]
    min_hours = input_time[3]; max_hours = input_time[9]
    min_mins = input_time[4]; max_mins = input_time[10]
    min_seconds = input_time[5]; max_seconds = input_time[11]
    #concatenate strings to form url
    url = 'https://tds.marine.rutgers.edu/thredds/ncss/roms/doppio/2017_da/his/History_Best?'\
    + 'var=u_eastward&var=v_northward&north=' + str(lat_max) + '&west=' + str(lon_min) + '&east='\
    + str(lon_max) + '&south=' + str(lat_min) + '&disableProjSubset=on&horizStride=1&time_start='\
    + str(min_year)+'-'+str(min_month)+'-'+str(min_day)+'T'+str(min_hours)+'%3A'+str(min_mins)+'%3A'+str(min_seconds)+'Z&time_end='\
    + str(max_year)+'-'+str(max_month)+'-'+str(max_day)+'T'+str(max_hours)+'%3A'+str(max_mins)+'%3A'+str(max_seconds)+'Z&timeStride=1&vertCoord='\
    + str(s_rho) + '&accept=netcdf'
    print(url)
    #url from website
    url = 'https://tds.marine.rutgers.edu/thredds/ncss/roms/doppio/2017_da/his/History_Best?var=u_eastward&var=v_northward&north=45.3302&west=-67.3578&east=-64.7255&south=44.2376&disableProjSubset=on&horizStride=1&time_start=2019-06-16T01%3A00%3A00Z&time_end=2019-06-19T12%3A00%3A00Z&timeStride=1&vertCoord=-0.9875&accept=netcdf'
    ds = netCDF4.Dataset(url) #make dataset
    print(ds)
    
    


def get_doppio_url(date):
    """ To fetch correct doppio url for a given time range 
    currently only for surface flow
    Written by Mark in 2021
    u is u_eastward
    v is v_northward
    time is time1
    lon is lon_rho, lat is lat_rho
    vertical levels are s_rho
    h is bathymetry at rho points
    -0.9875 is surface
    -0.0125 is bottom
    """
    url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/runs/History_RUN_2018-11T00:00:00Z'
    return url.replace('2018-11',date)





   
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

def bbox2ij(lon,lat,bbox):
    """Return indices for i,j that will completely cover the specified bounding box.
   
    i0,i1,j0,j1 = bbox2ij(lon,lat,bbox)
    
    lon,lat = 2D arrays that are the target of the subset
    bbox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]

    Example
    -------  
    >>> i0,i1,j0,j1 = bbox2ij(lon_rho,[-71, -63., 39., 46])
    >>> h_subset = nc.variables['h'][j0:j1,i0:i1]       
    """
    bbox=np.array(bbox)
   # np.set_printoptions(threshold=np.sys.maxsize)
    mypath=np.array([bbox[[0,1,1,0]],bbox[[2,2,3,3]]]).T
    #mypath_str = str(mypath)
    mypath_list=list(map(tuple,mypath))
    print(mypath)
    print(mypath_list)
    p = mpltpath.Path(mypath_list)
    points = np.vstack((lon.flatten(),lat.flatten())).T   
    n,m = np.shape(lon)
    inside = p.contains_points(points).reshape((n,m))
    ii,jj = np.meshgrid(range(m),range(n))
    return min(ii[inside]),max(ii[inside]),min(jj[inside]),max(jj[inside]) 



def rot2d(x, y, ang):
    '''rotate vectors by geometric angle
        
    This routine is part of Rob Hetland's OCTANT package:
    https://github.com/hetland/octant
    '''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


def shrink(a,b):
    """Return array shrunk to fit a specified shape by triming or averaging.
    
    a = shrink(array, shape)
    
    array is an numpy ndarray, and shape is a tuple (e.g., from
    array.shape). a is the input array shrunk such that its maximum
    dimensions are given by shape. If shape has more dimensions than
    array, the last dimensions of shape are fit.
    
    as, bs = shrink(a, b)
    
    If the second argument is also an array, both a and b are shrunk to
    the dimensions of each other. The input arrays must have the same
    number of dimensions, and the resulting arrays will have the same
    shape.
    
    This routine is part of Rob Hetland's OCTANT package:
    https://github.com/hetland/octant
    
    Example
    -------
    
    >>> shrink(rand(10, 10), (5, 9, 18)).shape
    (9, 10)
    >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
    [(5, 9, 10), (5, 9, 10)]   
       
    """

    if isinstance(b, np.ndarray):
        if not len(a.shape) == len(b.shape):
            raise Exception()
            'input arrays must have the same number of dimensions'
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)

    if isinstance(b, int):
        b = (b,)

    if len(a.shape) == 1:                # 1D array is a special case
        dim = b[-1]
        while a.shape[0] > dim:          # only shrink a
            if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                a = a[1:-1]
            else:                        # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])
    else:
        for dim_idx in range(-(len(a.shape)),0):
            dim = b[dim_idx]
            a = a.swapaxes(0,dim_idx)        # put working dim first
            while a.shape[0] > dim:          # only shrink a
                if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                    a = a[1:-1,:]
                if (a.shape[0] - dim) == 1:  # or average adjacent cells
                    a = 0.5*(a[1:,:] + a[:-1,:])
            a = a.swapaxes(0,dim_idx)        # swap working dim back

    return a


def plot_ROMS_velocity_field():
    """From Rich Signell's blog, plots ROMS velocity field for a single timepoint"""
    url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/runs/History_RUN_2018-05-15T00:00:00Z'
    nc = netCDF4.Dataset(url)
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    #bbox = [-71., -63.0, 41., 44.] #GoM
    bbox = [-67.35, -64.72, 44.23, 45.33] #BoF
    i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
    tvar = nc.variables['ocean_time']      # usual ROMS
    #tvar = nc.variables['time']     # USGS COAWST FMRC Aggregation
    h = nc.variables['h'][j0:j1, i0:i1]
    lon = lon_rho[j0:j1, i0:i1]
    lat = lat_rho[j0:j1, i0:i1]
    land_mask = 1 - nc.variables['mask_rho'][j0:j1, i0:i1]
    #start=datetime.datetime(2012,1,1,0,0)
    #start = datetime.datetime.utcnow()
    #tidx = netCDF4.date2index(start,tvar,select='nearest') # get nearest index to now
    tidx = -1
    #timestr = netCDF4.num2date(stats.tvar[tidx], stats.tvar.units).strftime('%b %d, %Y %H:%M') #BRING BACK
        
    zlev = -1  # last layer is surface layer in ROMS
    u = nc.variables['u'][tidx, zlev, j0:j1, i0:(i1-1)]
    v = nc.variables['v'][tidx, zlev, j0:(j1-1), i0:i1]
    
    lon_u = nc.variables['lon_u'][ j0:j1, i0:(i1-1)]
    lon_v = nc.variables['lon_v'][ j0:(j1-1), i0:i1]
    lat_u = nc.variables['lat_u'][ j0:j1, i0:(i1-1)]
    lat_v = nc.variables['lat_v'][ j0:(j1-1), i0:i1]
    
    lon=lon_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
    lat=lat_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
    mask = 1 - nc.variables['mask_rho'][(j0+1):(j1-1), (i0+1):(i1-1)]
    ang = nc.variables['angle'][(j0+1):(j1-1), (i0+1):(i1-1)]
    
    # average u,v to central rho points
    u = shrink(u, mask.shape)
    v = shrink(v, mask.shape)
    
    # rotate grid_oriented u,v to east/west u,v
    u, v = rot2d(u, v, ang)
    
    basemap = Basemap(projection='merc',llcrnrlat=44,urcrnrlat=46,llcrnrlon=-68,urcrnrlon=-64, lat_ts=30,resolution='i')
    fig1 = plt.figure(figsize=(10,8))
    ax = fig1.add_subplot(111)
    
    basemap.drawcoastlines()
    basemap.fillcontinents()
    basemap.drawcountries()
    basemap.drawstates()
    x_rho, y_rho = basemap(lon,lat)
    
    spd = np.sqrt(u*u + v*v)
    #h1 = basemap.pcolormesh(x_rho, y_rho, spd, vmin=0, vmax=1.0,shading='nearest') #add color
    nsub=2
    scale=0.03
    basemap.quiver(x_rho[::nsub,::nsub],y_rho[::nsub,::nsub],u[::nsub,::nsub],v[::nsub,::nsub],scale=1.0/scale, zorder=1e35, width=0.002)
    #basemap.colorbar(h1,location='right',pad='5%') #add colorbar
    title('COAWST Surface Current: ROMS Velocity Field') #BRING BACK
    plt.savefig('ROMS_velocity_field_BoF05152018.png')
    
def plot_ROMS_binned_velocity_field():
    "From uvmodel_function- IN DEVELOPMENT"
    #### THE ISSUE WITH THIS PROGRAM IS I HAVE MODEL OUTPUT FOR ONE DAY AND 
    ### I AM CALLING THE FUNCTION FOR SEVERAL DAYS... I NEED TO CHANGE GET_DOPPIO_URL
    date = '2018-05'
    url = get_doppio_url(date)
    nc = netCDF4.Dataset(url)
    depth=-1 #surface
    #depth=-15
    START_TIME= '2018-05-01 00:00:00' #CHANGE THIS
    END_TIME = '2018-05-15 00:00:00'
    START_TIME=dt.strptime(START_TIME, "%Y-%m-%d %H:%M:%S") 
    END_TIME=dt.strptime(END_TIME, "%Y-%m-%d %H:%M:%S") 
    # START_TIME = dt(2019, 5, 2, 0, 0)
    # END_TIME = dt(2019, 5, 17, 0, 0)
    timeperiod = (END_TIME)-(START_TIME)
    #startrecord=int(((timeperiod.total_seconds())/60/60))
    startrecord = 0
    
    
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    #bbox = [-71., -63.0, 41., 44.] #GoM
    bbox = [-67.5, -65, 44, 45.3] #BoF
    i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
    tvar = nc.variables['ocean_time']      # usual ROMS
    #tvar = nc.variables['time']     # USGS COAWST FMRC Aggregation
    h = nc.variables['h'][:]#[j0:j1, i0:i1]
    lon = lon_rho#[j0:j1, i0:i1]
    lat = lat_rho#[j0:j1, i0:i1]
    #land_mask = 1 - nc.variables['mask_rho']#[j0:j1, i0:i1]
    #start=datetime.datetime(2012,1,1,0,0)
    #start = datetime.datetime.utcnow()
    #tidx = netCDF4.date2index(start,tvar,select='nearest') # get nearest index to now
    tidx = -1
    #timestr = netCDF4.num2date(stats.tvar[tidx], stats.tvar.units).strftime('%b %d, %Y %H:%M') #BRING BACK
    
    zlev = -1  # last layer is surface layer in ROMS
    u = nc.variables['u_eastward'][tidx, zlev, j0:j1, i0:(i1-1)]
    v = nc.variables['v_northward'][tidx, zlev, j0:(j1-1), i0:i1]
    print(u)
    
    
    
    # lon_u = nc.variables['lon_u']#[ j0:j1, i0:(i1-1)]
    # lon_v = nc.variables['lon_v']#[ j0:(j1-1), i0:i1]
    # lat_u = nc.variables['lat_u']#[ j0:j1, i0:(i1-1)]
    # lat_v = nc.variables['lat_v']#[ j0:(j1-1), i0:i1]
    
    #lon=lon_rho[(j0+1):(j1-1), (i0+1):(i1-1)]
    #lat=lat_rho[(j0+1):(j1-1), (i0+1):(i1-1)]

    print("timeperiod= "+str(timeperiod))
    print("startrecord= "+str(startrecord))
    lon_bin_size= 0.05
    lat_bin_size= 0.05
    gbox=[-67.5, -65, 44, 45.3] #Bay of Fundy


    #siglay=nc.variables['siglay']

        
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

def plot_ROMS_depth():
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    bbox = [-71., -63.0, 41., 44.] #GoM
    i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,bbox)
    #tvar = nc.variables['ocean_time']      # usual ROMS
    tvar = nc.variables['time']     # USGS COAWST FMRC Aggregation
    h = nc.variables['h'][j0:j1, i0:i1]
    lon = lon_rho[j0:j1, i0:i1]
    lat = lat_rho[j0:j1, i0:i1]
    land_mask = 1 - nc.variables['mask_rho'][j0:j1, i0:i1]
    plt.pcolormesh(lon,lat,ma.masked_where(land_mask,-h),vmin=-350,vmax=0, shading='auto')
    plt.colorbar()
    plt.savefig('ROMSdepth.png')
    
    

###############MAIN PROGRAM############################
# time = '2019-12-30 12:00:00'
# if type(time)==str:
#     date_time=datetime.datetime.strptime(time,'%Y-%m-%d %H:%M:%S') # transform time format
# elif type(time)==datetime.datetime:
#     date_time=time
    
# #url_time=(date_time-datetime.timedelta(days=m)).strftime('%Y-%m-%d')#use this
# url_time=date_time.strftime('%Y-%m-%d')
# doppio_url=get_doppio_url(url_time) 
# nc = netCDF4.Dataset(doppio_url)



if __name__ == '__main__':
    getroms_field_byinfo(BoF_box,time_range,s_rho)
    
    
    
    