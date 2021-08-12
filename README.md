# Bay_of_Fundy
Assorted code from 2021 internship for simulated and observational drifter data extraction and analysis.


## Module Library

### Data Extraction/Visualizations
#### drifter_functions_2021.py
- Information: A collection of functions to extract observational drifter data from ERDDAP, plot the drifter on googlemaps and basemap.
- Example Output(s): 
  - Drifter Releases from 1988-2019 in Bay of Fundy

![Image of Drifter Releases from 1988-2019 in Bay of Fundy](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/drifterrelease_gmap_screenshot.png)


#### FVCOM_uvmodel.py
- Information: A collection of functions for visualizing FVCOM velocity fields, with tidal inputs and detided.
- Example Output(s):
  - FVCOM velocity field with tides (May 06-15 2016)

![Image of FVCOM velocity field with tides](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/FVCOM_velocity_field_with_tides_05062016-05152016_1m.png)

*Or you could detide the data*

  - FVCOM velocity field detided (May 01-15 2016)

![Image of FVCOM velocity field detided](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/FVCOM_velocity_field_detided_05012016-05152016_30m.png)


  - Subsurface currents not pictured


#### ROMS_uvmodel.py
- Information: Inspired by PICO github, uvmodel.py and Rich Signell's blog. Still working on binned and detided data processing.
  - ROMS COAWST Surface Current Gulf of Maine

![Image of COAWST Surface Current Gulf of Maine](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/ROMS_velocity_field_nocolor.png)

*Zoom into our area of interest*

  - ROMS COAWST Surface Current Bay of Fundy

![Image of COAWST Surface Current Bay of Fundy](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/ROMS_velocity_field_nocolor_BoF.png)

  - Detided plots coming soon


### Simulated Drifters
#### particle_tracking_routine
- Information: Using code from oceanparcels ipynb examples as posted on oceanparcels github. Solved with Runge-Kutta4 advection scheme.
  - DOPPIO Simulated Particles in BoF Gyre, May 2018

![Image of DOPPIO Simulated Particles](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/BoFParticleTrack1.png)

*What if we increase the number of particles? (sensitivity analysis coming soon)*

  - DOPPIO Simulated Particles in BoF Gyre, May 2018

![Image of DOPPIO Simulated Particles](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/BoFParticleTrack_05.01.2018_05.15.2018_10_2.png)

#### observed_vs_simulated
- Information: By advecting a particle in OceanParcels (RK4) and plotting the trajectory against an observed drifter released at the same lcoation over the same time period, generate this graph
  - DOPPIO Simulated particle vs observed drifter
 
![Image of DOPPIO Simulated Particles vs obs drifter](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/obs_vs_simdoppio_06_16_19_07_07_19.png)

Now, if we want to look at the separation distance (haversine) between those observed and simulated drifters (common model evalutation metric) for the first 72 hours:
![Image of DOPPIO Simulated Particles vs obs drifter sep distance](https://github.com/markeleone/Bay_of_Fundy/blob/main/plots/DOPPIO_sep_distance_over_time.png)

  
### Data Analysis/Statistics
- Coming soon: 
  - Sensitivity analysis
  - Separation distances (normalized)
  - Other Lagrangian Statistics
  - Forcing mechanisms
  - Additional plots
  - and more!
  
