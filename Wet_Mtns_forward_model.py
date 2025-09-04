# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from landlab.components import LinearDiffuser, FlowAccumulator, FastscapeEroder, DepressionFinderAndRouter
from landlab import RasterModelGrid, imshow_grid
from landlab.io import esri_ascii
from os import chdir

"""
This code was started on August 26, 2025.

The ultimate goal of this code is to use landlab to create a forward evolution model of the southern Wet Mountains, 
to visualize what this landscape looks like at equilibrium. We found evidence of drainage divide migration, as well as
stream capture. This model will hopefully elucidate any tectonic perturbations that could have initiated such landscape 
transience.

v 1. Create code that works
v 2. Download 30m DEM for Wets
v 3. Convert to ASCII (Turkey Creek from study site to confluence with Huerfano)
v 4.Test to see if works (use landscape evolution tutorials online landlab for imshow plots) "Quantifying river evolution with Landlab"
v 5. Look up uplift rates
v 6. Look up hillslope diffusivity
v 7. Make code for power law line to figure out m, n, and Ks (and refresh memory on Gallen 2013 and Gallen 2018)

"""

##instantiate grid

with open("basins_wets_30mbigclip.asc") as fp:
    mg = esri_ascii.load(fp, name = "topographic__elevation", at = "node")
z = mg.field_values('topographic__elevation')
mg.set_watershed_boundary_condition(z, -9999)
#mg.set_closed_boundaries_at_grid_edges(False, False, False, False) #set left and right closed, leave top and bottom open

#plt.figure(1)
#mg.imshow(mg.status_at_node, color_for_closed = "blue")

### ------- parameters / initial conditions
uplift_rate = 0.0001 #Marder et al, 2025 and Kainz et al 2023

n_sp = 4/3 # Whipple et al state that n = 1 for plucking and n = 5/3 for abrasion dominated erosion,
           # therefore n is probably somewhere in between (Whipple et al., 2000, Whipple & Tucker, 2002, Tucker & Whipple 2002)
theta = 0.27 #average of 0.31 for Turkey, 0.23 for St Charles found from 10m dem
m_sp = theta * n_sp 

K_sp = 1.54e-5 # found as average of Turkey Creek (5.8e-6) and St Charles basins (2.5e-5)

min_drainage = 1.25e5 # from Gallen et al., 2013 supplementary data
hillslope_diff = 0.005 #m2/year, from Pelletier et al., 2011; Langston et al., 2015; Struble et al., 2024

### ------ INSTANTIATE COMPONENTS
ld = LinearDiffuser(mg, linear_diffusivity = hillslope_diff)
fa = FlowAccumulator(mg, flow_director = 'D8', depression_finder = DepressionFinderAndRouter(mg))
fsc = FastscapeEroder(mg, K_sp = K_sp, m_sp = m_sp, n_sp = n_sp)

print('finished component instantiate')

active_time = 0
dt = 1000
total_time = 1000

while active_time < total_time:
    z[mg.core_nodes] += uplift_rate*dt
    print('finished uplift')
    ld.run_one_step(dt)
    print('finished diffusivity')
    fa.run_one_step()
    print('finished flow accumulate')
    fsc.run_one_step(dt)
    print('finished erosion')
    
    # if active_time % 5e6 == 0:
    #     plt.figure(3)
    #     mg.imshow("topographic__elevation", grid_units=("m", "m"), var_name="Elevation(m)")
    #     title_text = f"$K_{{sp}}$={K_sp}; $time$={total_time} yr; $dx$={30}m"
    #     plt.title(title_text)
    #     max_elev = np.max(z)
    #     print(f"Maximum elevation is {np.max(z)}")
        
    if active_time%100==0:
        print('finished time %s' % active_time)
    
    if active_time == 5e6:
         with open('WetsBig5Ma.asc', 'w') as f:
              esri_ascii.dump(mg, stream = f, at = 'node', name = 'topographic__elevation')
              
    if active_time == 10e6:
         with open('WetsBig10Ma.asc', 'w') as f:
              esri_ascii.dump(mg, stream = f, at = 'node', name = 'topographic__elevation')
              
    if active_time == 15e6:
         with open('WetsBig15Ma.asc', 'w') as f:
              esri_ascii.dump(mg, stream = f, at = 'node', name = 'topographic__elevation')
              
    if active_time == 20e6:
         with open('WetsBig20Ma.asc', 'w') as f:
              esri_ascii.dump(mg, stream = f, at = 'node', name = 'topographic__elevation')
    
    active_time += dt

print('finished model')          

