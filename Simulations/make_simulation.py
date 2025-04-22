
# Import modules
import numpy as np
from src.setup import New_Simulation

##############################
#       --- Inputs ---       #
# Please change this section #
##############################

# What do you want to call the simulation?
# This is just the name of the directory and does not affect the simulation.
# Can be set as None, and the script will be called something sensible based on the mesh.
simulation_name = 'P_prem_ani_source_0'

# What model do you want to use for the mesh?
# Either an in-built model or .bm file.
# In-built options include prem_iso, prem_ani, ak135, iasp91
model_name_for_mesh = 'prem_ani'

# What minimum period do you desire in seconds?
period_in_seconds = 20.

# Over what epicentral distance range do you want the database to cover?
# 0 to 180 is possible.
minimum_epicentral_distance = 0.
maximum_epicentral_distance = 90.

# Maximum depth of the database in km, this is deepest that sources can be.
# There isn't much point setting this deeper than 750 km as earthquakes don't happen
# that deep.
maximum_source_depth = 150.

# What order of polynomial to use for the database?
# At present can either be 2 or 4.
# Instaseis traditionally used 4, and this is most accurate, but 2 is almost as accurate
# and requires less storage. If you wish to use 2 then you will (at present) need to 
# edit instaseis.
polynomial_order = 4

# Source type for the simulation - must be PX, PZ, or both
# Use PX if you only want horizontal components, use PZ if you only want vertical, and
# use both if you want all 3 components.
source = 'both'

# How much seismogram do you want to calculate in seconds?
seismogram_length_in_seconds = 1500.

# Do you also want to run the mesher? True or False
# You might as well do this here, but you can do it yourself if you prefer
run_mesher = True

# Output writing buffer
# Reduce this if you run out of memory during the simulation.
buffer_size = 1000

########################################################
# Do not change this section unless you really want to #
########################################################

# Initialse simulation
New_Simulation(source=source, period=period_in_seconds, 
                seismogram_length=seismogram_length_in_seconds, 
                sim_name=simulation_name, mesh_file=model_name_for_mesh, 
                polynomial_order=polynomial_order, 
                min_distance=np.radians(minimum_epicentral_distance),
                max_distance=np.radians(maximum_epicentral_distance),
                max_depth_in_km=maximum_source_depth, 
                model_radius_in_km=None,
                source_magnitude=1e20,
                destination='.', draft_location='./Draft_Simulation/',
                run_mesher=run_mesher, buffer_size=buffer_size)
