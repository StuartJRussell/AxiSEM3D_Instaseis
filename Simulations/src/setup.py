
# Import modules
import os
import shutil
import logging
import subprocess
import ruamel.yaml
import numpy as np

# Load YAML and set some formatting settings
yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml.default_flow_style = None
yaml.indent(mapping=4, sequence=4, offset=4)

# Logger and configuration
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class New_Simulation:

    def __init__(self, source, period, seismogram_length, 
                        sim_name=None, mesh_file='prem_iso', polynomial_order=4,
                        min_distance=0., max_distance=np.pi,
                        max_depth_in_km=750., model_radius_in_km=None,
                        source_magnitude = 1e20,
                        destination='.', draft_location='./Draft_Simulation/',
                        run_mesher=False, buffer_size=1000):
        
        # Ensure either PX, PZ or both
        if source not in ['PX', 'PZ', 'both']:
            raise ValueError('Source must be PX, PZ or both')
        
        # Ensure polynomial order is either 2 or 4
        if polynomial_order not in [2, 4]:
            raise ValueError('Polynomial order must be 2 or 4')
        
        # Add inputs as attributes
        self.mesh_file = mesh_file
        self.mesh_file_location = mesh_file
        self.mesh_name = self.mesh_file
        self.source = source
        self.source_name = source if source != 'both' else 'P'
        self.period = period
        self.seismogram_length = int(seismogram_length)
        self.period_str = f'{period:.0f}s'
        self.polynomial_order = polynomial_order
        self.min_distance = np.max([0., min_distance])
        self.max_distance = np.min([max_distance, np.pi])
        self.max_depth = max_depth_in_km * 1e3
        self.model_radius = model_radius_in_km
        self.source_magnitude = source_magnitude
        self.draft_location = draft_location
        self.run_mesher = run_mesher
        self.buffer_size = buffer_size
        
        # If mesh is .bm file then get details
        # At the moment just the mesh name
        if '.bm' in self.mesh_file:
            self.get_mesh_details(self.mesh_file)
        
        # Decide simulation name
        if sim_name is None:
            sim_name = '_'.join([self.source_name, self.mesh_name, self.period_str])
        self.sim_name = sim_name
        self.final_destination = os.path.abspath(os.path.join(destination, sim_name))
        
        # Copy files over to the new location
        directory_made = self.make_destination()
        
        # Paths to input files
        self.inparam_advanced = os.path.join(
                            self.final_destination, 'input', 'inparam.advanced.yaml')
        self.inparam_mesh = os.path.join(
                            self.final_destination, 'input', 'inparam.mesh.yaml')
        self.inparam_model = os.path.join(
                            self.final_destination, 'input', 'inparam.model.yaml')
        self.inparam_nr = os.path.join(
                            self.final_destination, 'input', 'inparam.nr.yaml')
        self.inparam_output = os.path.join(
                            self.final_destination, 'input', 'inparam.output.yaml')
        self.inparam_source = os.path.join(
                            self.final_destination, 'input', 'inparam.source.yaml')

        # Set the mesh inputs file
        logger.info('  Dealing with the mesh')
        self.write_mesh_inputs()
        
        # Run the mesher if required
        if run_mesher:
            self.create_mesh()
        else:
            logger.warning('    I am not running the mesher, you need to do this yourself')
        
        # Edit the simulation input files
        self.edit_input_files()
        
        logger.info(f'Simulation created')
        


    def make_destination(self):
        """
        This function copies the draft simulation to a new directory.
        """
        
        # If destination already exists then ask for advice
        if os.path.exists(self.final_destination):
            answer = input('Destination already exists, shall I overwrite it? y/n: ')
            if answer != 'y':
                logger.warning('Okay, not overwriting, quitting instead')
                exit()
            else:
                shutil.rmtree(self.final_destination)
        
        # Copy files
        shutil.copytree(self.draft_location, self.final_destination)
        if '.bm' in self.mesh_file:
            shutil.copy(self.mesh_file, self.final_destination + '/input/')
        logger.info('Simulation directory created:')
        logger.info(f'    {self.final_destination}')
        
        return True



    def write_mesh_inputs(self):
        """
        This function edits inparam.mesh.yaml.
        """
        
        logger.info('    Creating mesh inputs')
        
        # Open file
        filename = self.inparam_mesh
        with open(filename) as f:
            parameters = yaml.load(f)
        
        # Change the parameters
        parameters['basic']['model'] = f'{self.mesh_file_location}'
        parameters['basic']['period'] = self.period
        
        # Write file
        with open(filename, "w") as f:
            yaml.dump(parameters, f)
        
        # Create exodus file name
        if not self.run_mesher:
            self.exodus_file= '_'.join([parameters['mesh_type'], 
                                         self.mesh_name,
                                         self.period_str]) + '.e'



    def edit_input_files(self):
        """
        Function just edits the input files in turn.
        """
        
        logger.info('  Setting the input files')
    
        # Set the model inputs file
        self.write_model_inputs()
        
        # Set the nr inputs file
        self.write_nr_inputs()
        
        # Set the output inputs file
        self.write_output_inputs()
        
        # Set the output inputs file
        self.write_source_inputs()
        
        logger.info('  Finished setting the input files')



    def write_model_inputs(self):
        """
        This function edits inparam.model.yaml.
        """
        
        logger.info('    Creating model inputs')
        
        # Open file
        filename = self.inparam_model
        with open(filename) as f:
            parameters = yaml.load(f)
        
        # Set the 1D model
        parameters['model1D']['exodus_mesh'] = self.exodus_file
        
        # Geodesy settings
        parameters['geodesy']['lat_lon_north_pole_mesh'] = 'SOURCE'
        parameters['geodesy']['flattening_on_surface'] = 'SPHERE'
        
        # Ensure no 3D models
        parameters['list_of_3D_models'] = []

        # Write file
        with open(filename, "w") as f:
            yaml.dump(parameters, f)



    def write_nr_inputs(self):
        """
        This function edits inparam.nr.yaml.
        """
        
        logger.info('    Creating nr inputs')
        
        # Open file
        filename = self.inparam_nr
        with open(filename) as f:
            parameters = yaml.load(f)
        
        # Set Nr values
        parameters['type_Nr'] = 'CONSTANT'
        parameters['constant'] = 1 if self.source == 'PZ' else 3

        # Write file
        with open(filename, "w") as f:
            yaml.dump(parameters, f)



    def write_output_inputs(self):
        """
        This function edits inparam.output.yaml.
        """
        
        logger.info('    Creating output inputs')
        
        # Take or guess the radius
        if hasattr(self, 'mesh_radius'):
            mesh_radius = self.mesh_radius
            logger.info(f'      Taking mesh radius of {mesh_radius} km from mesh file')
        elif self.model_radius is not None:
            mesh_radius = self.model_radius * 1e3
            logger.info(f'      Taking mesh radius of {mesh_radius} km as inputted')
        else:
            mesh_radius = 6371e3
            logger.info(f'      Assuming mesh radius of {mesh_radius} km')
        
        # Open file
        filename = self.inparam_output
        with open(filename) as f:
            parameters = yaml.load(f)
        oparameters = parameters.copy()
        
        # Set spatial ranges of the output
        parameters['list_of_element_groups'][0]['mantle'][
                                'elements']['horizontal_range'] = [
                                            float(self.min_distance),
                                            float(self.max_distance)]
        parameters['list_of_element_groups'][0]['mantle'][
                                'elements']['vertical_range'] = [
                                            int(mesh_radius - self.max_depth),
                                            int(mesh_radius)]
        
        # Inplance parameters
        parameters['list_of_element_groups'][0]['mantle'][
                        'inplane']['edge_dimension'] = 'BOTH'
        if self.polynomial_order == 4:
            parameters['list_of_element_groups'][0]['mantle'][
                            'inplane']['GLL_points_one_edge'] = [0, 1, 2, 3, 4]
        elif self.polynomial_order == 2:
            parameters['list_of_element_groups'][0]['mantle'][
                            'inplane']['GLL_points_one_edge'] = [0, 2, 4]
        
        # Azimuthal parameters
        parameters['list_of_element_groups'][0]['mantle'][
                        'azimuthal']['phi_list'] = []
        parameters['list_of_element_groups'][0]['mantle'][
                        'azimuthal']['lat_lon_list'] = []
        parameters['list_of_element_groups'][0]['mantle'][
                        'azimuthal']['na_space'] = 1
        
        # Wavefield parameters
        parameters['list_of_element_groups'][0]['mantle'][
                        'wavefields']['coordinate_frame'] = 'spz'
        parameters['list_of_element_groups'][0]['mantle'][
                        'wavefields']['medium'] = 'SOLID'
        parameters['list_of_element_groups'][0]['mantle'][
                        'wavefields']['channels'] = ['U']
        
        # Temporal
        parameters['list_of_element_groups'][0]['mantle'][
                        'temporal']['sampling_period'] = 'DTx40'
        parameters['list_of_element_groups'][0]['mantle'][
                        'temporal']['time_window'] = 'FULL'
        
        # Buffer
        parameters['list_of_element_groups'][0]['mantle'][
                        'file_options']['buffer_size'] = self.buffer_size

        # Write file
        with open(filename, "w") as f:
            yaml.dump(parameters, f)



    def write_source_inputs(self):
        """
        This function edits inparam.source.yaml.
        """
        
        logger.info('    Creating source inputs')
        
        # Open file
        filename = self.inparam_source
        with open(filename) as f:
            parameters = yaml.load(f)
        
        # Set the length of the seismogram
        parameters['time_axis']['record_length'] = self.seismogram_length
        
        # Set the source location
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE'][
                                    'location']['latitude_longitude'] = [90e0, 0e0]
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE'][
                                    'location']['depth'] = 0e0
        
        # Set the mechanism
        force_vectors = {'PZ':[1e0, 0e0, 0e0],
                         'PX':[0e0, 1e0, 0e0],
                         'both':[1e0, 1e0, 0e0]}
        force_vector = force_vectors[self.source]
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE']['mechanism'][
                                        'type'] = 'FORCE_VECTOR'
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE']['mechanism'][
                                        'data'] = force_vector
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE']['mechanism'][
                                        'unit'] = float(self.source_magnitude /
                                         np.linalg.norm(force_vector))
        
        # Source time function
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE'][
                    'source_time_function']['class_name'] = 'GaussianSTF'
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE'][
                    'source_time_function']['half_duration'] = self.period
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE'][
                    'source_time_function']['time_shift'] = 0e0
        parameters['list_of_sources'][0]['THE_ONLY_SOURCE'][
                    'source_time_function']['use_derivative_integral'] = 'GAUSSIAN'
        
        # Write file
        with open(filename, "w") as f:
            yaml.dump(parameters, f)



    def get_mesh_details(self, mesh_file):
        """
        Reads the .bm and gets the mesh details.
        """
        
        with open(mesh_file, 'r') as f:
            lines = f.readlines()
            name_line = [x for x in lines if 'NAME' in x][0]
            
            # Set name from the file
            self.mesh_name = name_line.strip().split()[1]
            
            # Location for salvusmeshlite
            self.mesh_file_location = './input/' + self.mesh_file.split('/')[-1]



    def create_mesh(self):
        """
        This function simply runs the makemeshgrid.sh script.
        The user can of course do this themselves if they prefer.
        """
        
        # Get where we currently are
        start_dir = os.getcwd()
        
        # Navigate into the simulation directory
        os.chdir(self.final_destination)
        
        # Run makemeshgrid
        logger.info('    Running the mesher')
        out = subprocess.run(["bash", "-l", "makemeshgrid.sh"], capture_output=True,
                                text=True)
        if out.returncode != 0:
            raise Exception('Mesher has failed')
        self.exodus_file = out.stdout.split('\n')[-2].split("'")[1]
        self.mesh_radius = float([
                x for x in out.stdout.split('\n') if 'radius' in x][0].split()[-1])
        
        # Navigate back
        os.chdir(start_dir)
        logger.info('    Mesher has run')





