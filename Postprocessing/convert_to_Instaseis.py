
# Import modules
import sys
import time
import logging
import numpy as np

# Import class
from src.core import AxiSEM3D_Output

# Logger and configuration
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

##############################
#       --- Inputs ---       #
# Please change this section #
##############################

# Which simulation
sim_dir = '../Simulations/'
sim_names = ['P_prem_ani_20s']

# Path to output database
database_dir = '../Databases/'
database_name = 'AxiSEM3D_prem_ani_20s'
path = f'{database_dir}{database_name}/'

########################################################
# Do not change this section unless you really want to #
########################################################

t0 = time.time()
processes = int(sys.argv[1])
for sim_name in sim_names:

    # Initialise object
    p = AxiSEM3D_Output(f'{sim_dir}/{sim_name}')
    
    # Create database
    p.create_instaseis_database(path, processes=processes)

t1 = time.time()
logger.info(f'Database creation took {t1 - t0} seconds with {processes} processes')
sys.exit()
