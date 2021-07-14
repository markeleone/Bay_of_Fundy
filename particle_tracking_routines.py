# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 11:36:34 2021

@author: markl
Using code from oceanparcels ipynb examples as posted on oceanparcels github
Used for BoF manuscript project by Mark in 2021
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
import netCDF4
import os
import pathlib as path
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj") #needed to run basemap
from mpl_toolkits.basemap import Basemap
from scipy import stats
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile

"""
#Example from oceanparcels ipynb, runs nicely except output file saving
filenames = {'U': "parcels_examples/GlobCurrent_example_data/20*.nc",
             'V': "parcels_examples/GlobCurrent_example_data/20*.nc"}
variables = {'U': 'eastward_eulerian_current_velocity',
             'V': 'northward_eulerian_current_velocity'}
dimensions = {'lat': 'lat',
              'lon': 'lon',
              'time': 'time'}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
pset = ParticleSet.from_line(fieldset=fieldset, pclass=JITParticle,
                             size=5,            # releasing 5 particles
                             start=(28, -33),   # releasing on a line: the start longitude and latitude
                             finish=(30, -33))  # releasing on a line: the end longitude and latitude
output_file = pset.ParticleFile(name="parcels_examples/GlobCurrentParticles.nc", outputdt=timedelta(hours=6))
pset.execute(AdvectionRK4,
             runtime=timedelta(days=10),
             dt=timedelta(minutes=5),
             output_file=output_file)
output_file.export()
plotTrajectoriesFile('parcels_examples/GlobCurrentParticles.nc',
                     tracerfile='parcels_examples/GlobCurrent_example_data/20020101000000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc',
                     tracerlon='lon',
                     tracerlat='lat',
                     tracerfield='eastward_eulerian_current_velocity');
"""


#url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/runs/History_RUN_2018-11T00:00:00Z'
#nc = netCDF4.Dataset(url)

######MAIN CODE FOR DOPPIO SIMULATIONS#################
filenames = {'U': 'DOPPIOOutput5.01.2018-5.15.2018.nc',
             'V': 'DOPPIOOutput5.01.2018-5.15.2018.nc'}
variables = {'U': 'u_eastward',
             'V': 'v_northward'}
dimensions = {'lat':'lat_rho',
              'lon':'lon_rho',
              'time':'time1'}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
pset = ParticleSet.from_line(fieldset=fieldset, pclass=JITParticle,
                             size = 5, #releasing x particles
                             start=(-66.8,45), #releasing on a line: the start longitude and latitude
                             finish=(-66.2,44.6)) # releasing on a line: the end longitude and latitude
output_file = pset.ParticleFile(name='BoFParticles.nc', outputdt=timedelta(hours=6))
pset.execute(AdvectionRK4, #4th order RK
             runtime=timedelta(days=14), #run for 10 days
             dt=timedelta(minutes=5),
             output_file = output_file)

#print('BoFParticles.nc')

output_file.export()

plt = plotTrajectoriesFile('BoFParticles.nc')
                     # tracerfile='his_History_Best.nc',
                     # tracerlon='lon_rho',
                     # tracerlat='lat_rho',
                     # tracerfield='ubar_eastward')"""

