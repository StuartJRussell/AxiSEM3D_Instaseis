
# Import modules
import os
import yaml
import logging
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
import multiprocessing
from scipy.special import erf
from multiprocessing import Process, Queue

# Import local files
from src import mesh
from src import tools

# Logger and configuration
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class AxiSEM3D_Output:
    """
    Class to handle and convert the element output of reciprocal AxiSEM3D simulations
    such that they can be used with Instaseis.

    At present some parts have been cannibalised from `axikernels` by Adrian Mag, these
    are indicated in the docstrings. These have largely been simplified due to 
    redundancies associated with only dealing with 1D models and force sources.
    
    Several functions in this package reference equations in Nissen-Meyer et al. 2007,
    the full reference is:
        Nissen‐Meyer, T., Fournier, A., & Dahlen, F. A. (2007). A two‐dimensional
        spectral‐element method for computing spherical‐earth seismograms–I. 
        Moment‐tensor source. Geophysical Journal International, 168(3), 1067-1092.
        doi: 10.1111/j.1365-246X.2006.03121.x
    
    STILL TO DO:
    - Check amplitudes of ERF sources
    - Write merged database rather than separate?
    """
    
    def __init__(self, simulation_path, input_dir='input', 
                                element_dir='output/elements',
                                output_all_mesh_parameters=False):
        
        # Path to simulation
        full_simulation_path = os.path.abspath(simulation_path)
        if full_simulation_path[-1] != '/':
            full_simulation_path += '/'
        
        # Set some attributes based on the inputs
        self.path_to_simulation = full_simulation_path
        self.path_to_inputs = ''.join([self.path_to_simulation, input_dir, '/'])
        self.path_to_elements = ''.join([self.path_to_simulation, element_dir, '/'])
        self.output_all_mesh_parameters = output_all_mesh_parameters
        
        # Read the input files
        self._load_input_files()
        
        # Create the metadata for the elements
        self._load_elements_info()
        
        # Set chunk size to be a little over the length of the number of GLL points
        self.chunksize = len(self.element_group_info['inplane']['sem_inds']) + 2



    def _load_input_files(self):
        """
        This function reads in the inparam files associated with the AxiSEM3D run.
        It also then checks whether these have the corrcet settings.
        """
        
        # Define the paths to the original inparam files
        self.inparam_advanced_file = (
                    self.path_to_inputs + '/inparam.advanced.yaml')
        self.inparam_model_file = (
                    self.path_to_inputs + '/inparam.model.yaml')
        self.inparam_nr_file = (
                    self.path_to_inputs + '/inparam.nr.yaml')
        self.inparam_output_file = (
                    self.path_to_inputs + '/inparam.output.yaml')
        self.inparam_source_file = (
                    self.path_to_inputs + '/inparam.source.yaml')
        
        # Read the inparam files
        with open(self.inparam_advanced_file, 'r') as f:
            self.inparam_advanced=yaml.load(f, Loader=yaml.FullLoader)
        with open(self.inparam_model_file, 'r') as f:
            self.inparam_model=yaml.load(f, Loader=yaml.FullLoader)
        with open(self.inparam_nr_file, 'r') as f:
            self.inparam_nr=yaml.load(f, Loader=yaml.FullLoader)
        with open(self.inparam_output_file, 'r') as f:
            self.inparam_output=yaml.load(f, Loader=yaml.FullLoader)
        with open(self.inparam_source_file, 'r') as f:
            self.inparam_source=yaml.load(f, Loader=yaml.FullLoader)
        
        # Add mesh filename explicitly
        self.mesh_file = (self.path_to_inputs +
                            self.inparam_model['model1D']['exodus_mesh'])
        
        # Check that inputs are compatible with Instaseis
        self._check_inputs()



    def _check_inputs(self):
        """
        This function checks the inputs to make sure they are correct and compatible.
        """
        
        # Check source
        self._check_source()
        
        # Check element output
        self._check_output()
        
        # Check Nr
        self._check_Nr()



    def _check_source(self):
        """
        Checks the source inputs.
        """
        
        # Check number of sources
        if len(self.inparam_source['list_of_sources']) != 1:
            raise ValueError('More than one source')
        
        # Assign source information separately
        sname = list(self.inparam_source['list_of_sources'][0].keys())[0]
        self.source_info = self.inparam_source['list_of_sources'][0][sname]
        
        # Check location is at north pole
        source_latitude = float(self.source_info['location']['latitude_longitude'][0])
        source_longitude = float(self.source_info['location']['latitude_longitude'][1])
        if source_latitude != 90.:
            raise ValueError('Source location must be at the North Pole')
        
        # Check that the source is at the surface
        if float(self.source_info['location']['depth']) != 0.:
            raise ValueError('Source must have colatitude of 0')
        
        # Check mechanism is force
        if self.source_info['mechanism']['type'] != 'FORCE_VECTOR':
            raise ValueError('Source type must be FORCE_VECTOR')
        
        # Source vector
        source_vector = np.array(
                    [float(x) for x in self.source_info['mechanism']['data']])
        
        # Classify source type
        if len(source_vector[source_vector != 0.]) == 1:
            self.source_type = np.array(['PZ', 'PX0', 'PX1'])[source_vector != 0.][0]
        elif (source_vector[0] != 0. and source_vector[1] != 0. and 
              source_vector[2] == 0.):
             self.source_type = 'P'
        else:
            raise ValueError('Source vector not PX, PX or both.')

        # The magnitude is actually the magnitude of a single component
        # This is a left over peculiarity from AxiSEM calculating PX and PZ
        # separately.
        if self.source_type == 'P' and source_vector[0] != source_vector[1]:
            raise ValueError('For source type both, PX and PZ must be equal')
        if self.source_type == 'P':
            self.source_magnitude = (source_vector[0]
                        * float(self.source_info['mechanism']['unit']))
        else:
            self.source_magnitude = (np.linalg.norm(source_vector)
                        * float(self.source_info['mechanism']['unit']))
        
        # Ensure STF is Gaussian
        if self.source_info['source_time_function']['class_name'] != 'GaussianSTF':
            raise ValueError('Source time function must be GaussianSTF')



    def _check_output(self):
        """
        Checks the output inputs.
        """
        
        # Check that there is only one element group
        if len(self.inparam_output['list_of_element_groups']) > 1:
            raise ValueError('There is more than 1 element group')
        
        # Set the name of the group
        self.element_group_name = list(
                        self.inparam_output['list_of_element_groups'][0].keys())[0]
        
        # Set the information to an attribute that can be edited
        self.element_group_info = self.inparam_output[
                        'list_of_element_groups'][0][self.element_group_name].copy()
        
        # Check that element centre is outputted
        if 2 not in self.element_group_info['inplane']['GLL_points_one_edge']:
            raise ValueError('Element centre must be outputted')
        
        # What is the polynomial order based on the output
        # This is regardless of what was used for the simulation
        edge_output = self.element_group_info['inplane']['GLL_points_one_edge']
        npol = len(edge_output) - 1
        self.element_group_info['inplane']['npol'] = npol
        
        # Indices for mapping
        # Instaseis/AxiSEM assumes that points are in a certain spatial arrangement
        if edge_output == [2]:
            self.element_group_info['inplane']['mp_ind'] = 0
            self.element_group_info['inplane']['fem_inds'] = np.array([0])
            self.element_group_info['inplane']['sem_inds'] = np.array([0])
        elif edge_output == [0, 2, 4]:
            self.element_group_info['inplane']['mp_ind'] = 4
            self.element_group_info['inplane']['fem_inds'] = np.array([0, 6, 8, 2])
            self.element_group_info['inplane']['sem_inds'] = np.array([0, 3, 6, 1, 4,
                                                                       7, 2, 5, 8])
        elif edge_output == [0, 1, 2, 3, 4]:
            self.element_group_info['inplane']['mp_ind'] = 12
            self.element_group_info['inplane']['fem_inds'] = np.array([0, 20, 24, 4])
            self.element_group_info['inplane']['sem_inds'] = np.array([
                                0, 5, 10, 15, 20, 1, 6, 11, 16, 21, 2, 7, 12, 
                                17, 22, 3, 8, 13, 18, 23, 4, 9, 14, 19, 24])
        
        # Check that the medium is solid - needed for the mantle...
        if self.element_group_info['wavefields']['medium'] != 'SOLID':
            raise ValueError('Medium must be SOLID...')
        
        # Check the coordinate frame is spz
        if self.element_group_info['wavefields']['coordinate_frame'] != 'spz':
            raise ValueError('Coordinate_frame must be set to spz')
        
        # Check that displacement is recorded
        if 'U' not in self.element_group_info['wavefields']['channels']:
            raise ValueError('U channels must be recorded')
        
        # Check that sampling rate is relative to DT
        if 'DTx' not in self.element_group_info['temporal']['sampling_period']:
            raise ValueError('Database DT must be defined relative to simulation DT')



    def _check_Nr(self):
        """
        Checks that Nr is sufficient.
        Higher values of Nr are unnecessary but permissable.
        """
        
        # Check Nr is constant
        if self.inparam_nr['type_Nr'] != 'CONSTANT':
            raise ValueError('Nr is not constant')
        
        # Check that the value is high enough
        self.Nr = self.inparam_nr['constant']
        if self.source_type == 'PZ' and self.Nr < 1:
            raise ValueError('Nr must be at least 1 for source PZ')
        elif self.source_type in ['P', 'PX0', 'PX1'] and self.Nr < 3:
            raise ValueError('Nr must be at least 3 for sources involving PX')



    def _load_elements_info(self):
        """
        Adds metadata about the element group.
        Courtesy of Adrian Mag.
        """
        
        # Get and assign the metadata
        metadata = self._read_element_metadata(self.element_group_name)
        self.element_group_info['metadata'] = {
                    'na_grid': metadata[0],
                    'data_time': metadata[1],
                    'list_element_na': metadata[2],
                    'list_element_coords': metadata[3],
                    'list_element': metadata[4],
                    'files': metadata[5],
                    'elements_index_limits': metadata[6],
                    'detailed_channels': metadata[7]}
        
        # Set some variables that come from the metadata
        time_axis = self.element_group_info['metadata']['files'][0].data_time.data
        self.len_time = len(time_axis)
        self.dt = time_axis[1] - time_axis[0]
        self.start_time = time_axis[0]
        self.npoints = len(self._get_unique_points(return_all=False))
        self.nelements = self.element_group_info[
                                    'metadata']['elements_index_limits'][-1]



    def _read_element_metadata(self, element_group):
        """
        Collates the contents of the element output files.
        Courtesy of Adrian Mag.
        """
        
        # Filenames of the outputs
        path_to_element_group = os.path.join(self.path_to_elements,
                                             element_group)
        nc_fnames = [os.path.join(path_to_element_group, f)
                         for f in os.listdir(path_to_element_group) if
                            'axisem3d_synthetics.nc' in f]
        
        # Open files
        nc_files = []
        for nc_fname in nc_fnames:
            nc_files.append(xr.open_dataset(nc_fname))

        # Read Na grid (all azimuthal dimensions)
        # Only one value of Na for these simulations
        na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)

        # Read time axis
        data_time = nc_files[0].data_vars['data_time'].values

        # Variables must be concatenated over the datasets minus the data itself
        # Define some empty lists of xarray.DataArray objects
        xda_list_element_na = []
        xda_list_element_coords = []
        xda_list_element = []

        # Get channel names as letters
        coordinate_frame = self.element_group_info['wavefields']['coordinate_frame']
        detailed_channels = [element.replace('1', coordinate_frame[0])
                                    .replace('2', coordinate_frame[1])
                                    .replace('3', coordinate_frame[2])
                             for element in [str_byte.decode('utf-8') for
                                str_byte in nc_files[0].list_channel.data]]

        # Loop over nc files and append DataArrays
        index_limit = 0
        elements_index_limits = [0]
        for nc_file in nc_files:
            index_limit += nc_file.sizes['dim_element']
            elements_index_limits.append(index_limit)
            xda_list_element_na.append(nc_file.data_vars['list_element_na'])
            xda_list_element_coords.append(nc_file.data_vars['list_element_coords'])
            xda_list_element.append(
                        nc_file.data_vars['list_element__NaG=%d' % self.Nr])

        # Concatenate DataArrays
        xda_list_element_na = xr.concat(xda_list_element_na, dim='dim_element')
        xda_list_element_coords = xr.concat(xda_list_element_coords,
                                                dim='dim_element')
        xda_list_element = xr.concat(xda_list_element,
                                                   dim='dim_element__NaG=%d' % self.Nr)

        # Convert to numpy arrays
        list_element_na = xda_list_element_na.values.astype(int)
        list_element_coords = xda_list_element_coords.values.astype(np.float32)
        list_element = xda_list_element.values.astype(int)

        # return: Here we return the files only because in this format they are
        # not being loaded into RAM Since these files are huge we prefer to
        # load into RAM only the file where the data that we want is located
        # and then close the file.
        return na_grid, data_time, list_element_na, \
                list_element_coords, list_element, \
                nc_files, elements_index_limits, \
                detailed_channels
    
    
    
    def _get_all_coordinates(self):
        """
        Function gets all coordinates from the output files.
        # These are rounded to negate the effect of very small differences.
        """

        # All coordinates
        # Round to negate small differences
        coordinates = np.round(self.element_group_info['metadata'][
                        'list_element_coords'][:, :, :].reshape(-1, 2), 3)
        
        return coordinates
        


    def _get_unique_points(self, return_all):
        """
        Function returns the unique coordinates from within the element group and
        sorts according to a measure.
        """
        
        # Get all coordinates
        coordinates = self._get_all_coordinates()
        
        # Add theta to the array so that np.unique sorts by this
        theta = np.arctan2(coordinates[:,0], coordinates[:,1])
        r = np.sqrt(coordinates[:,0]**2. + coordinates[:,1]**2.)
        sort = np.round((16. * theta / np.pi) + 0.5) + ((np.max(r) - r) / np.max(r))
        coordinates = np.vstack((sort, coordinates.T)).T
        
        # Get the unique coordinates
        coordinates, unique_indices, unique_inverse = np.unique(coordinates, axis=0,
                        return_index=True, return_inverse=True)
        
        # Remove sorting array from coordinates before returning
        coordinates = coordinates[:,1:]
        
        # Return as requested
        if return_all:
            return coordinates, unique_indices, unique_inverse
        else:
            return coordinates



    def _get_element_mapping(self, what):
        """
        When removing repeated GLL points it becomes necessary to define how GLL points
        and elements map to each other. This function creates that mapping.
        """
        
        # Mids, FEM or SEM
        if what not in ['Mids', 'FEM', 'SEM']:
            raise ValueError('Mappping request must be for Mids, FEM or SEM')
        
        # Some parameters about the mesh
        # These all already exist, it just neatens the code
        npol = self.element_group_info['inplane']['npol']
        mp_ind = self.element_group_info['inplane']['mp_ind']
        fem_inds = self.element_group_info['inplane']['fem_inds']
        sem_inds = self.element_group_info['inplane']['sem_inds']
        nelements = self.element_group_info['metadata']['elements_index_limits'][-1]
        nGLL = np.prod(self.element_group_info[
                        'metadata']['list_element_coords'].shape[:2])

        # Take only unique coordinates
        coordinates, unique_indices, unique_inverse = self._get_unique_points(
                                                        return_all=True)
        
        # Midpoint mapping
        if what == 'Mids':
            mapping = unique_inverse[mp_ind::len(sem_inds)]
        
        # FEM mapping
        elif what == 'FEM':
            mapping = unique_inverse[np.array(
                                            [(fem_inds + i * ((npol + 1) ** 2)) 
                                                for i in range(nelements)])]
        
        # SEM mapping
        elif what == 'SEM':
            mapping = unique_inverse[np.array(
                                    [(sem_inds + i * ((npol + 1) ** 2)).reshape(
                                                            npol + 1, npol + 1)
                                                    for i in range(nelements)])]
        
        return mapping



    def create_instaseis_database(self, path, processes):
        """
        This function just calls the _create_netCDF_database with the appropriate
        source types.
        """
        
        if self.source_type == 'PZ':
            self._create_netCDF_database(path, 'PZ', processes)
        elif self.source_type == 'PX0':
            self._create_netCDF_database(path, 'PX0', processes)
        elif self.source_type == 'PX1':
            self._create_netCDF_database(path, 'PX1', processes)
        elif self.source_type == 'P':
            self._create_netCDF_database(path, 'PZ', processes)
            self._create_netCDF_database(path, 'PX0', processes)
    
    
    
    def _create_netCDF_database(self, path, source_type, processes):
        """
        This function creates and writes the netCDF database file.
        """
        
        # Ensure path is absolute
        path = os.path.abspath(path)
        
        # Create folder structure
        db_path = self._create_instaseis_structure(path, source_type)

        # Create netCDF file
        filename = db_path + '/ordered_output.nc4'
        
        # Open a netCDF dataset
        logger.info(f'  Creating netCDF file for source type {source_type}')
        nc_out = nc.Dataset(filename, 'w', format='NETCDF4')
        
        # Create the dimensions, only two required
        nc_out.createDimension('gllpoints_all', self.npoints)
        nc_out.createDimension('snapshots', self.len_time)
        
        # Create snapshots group
        self._create_output_snapshots(nc_out, source_type, processes)
        
        # Create Mesh group
        self._create_output_mesh(nc_out)
        
        # Create the attributes for the file
        attributes = self._create_database_attributes(source_type)
        for attr in attributes:
            setattr(nc_out, attr, attributes[attr])
        
        # Save and close file
        nc_out.close()
        logger.info(f'{source_type} database written: {filename}')



    def _create_instaseis_structure(self, path, source_type):
        """
        Simply creates the output directory structure where the database will be stored.
        """
        
        # Paths to database
        db_path = os.path.join(path, source_type[:2], 'Data')
        
        # Make output folder if it doesn't already exist
        if not os.path.exists(db_path):
            os.makedirs(db_path)
        
        return db_path



    def _create_output_snapshots(self, nc_out, source_type, processes):
        """
        Creates the `Snapshots` group which contains the displacement and STF data.
        """
        
        # Create the group
        nc_out.createGroup('Snapshots')
    
        # Add nstrain attribute to Snapshots
        setattr(nc_out.groups['Snapshots'], 'nstrain', self.len_time)
        
        # Create disp_s variable
        channels = self.element_group_info['metadata']['detailed_channels']
        disp_s = nc_out.groups['Snapshots'].createVariable("disp_s", np.float32,
                                            ("snapshots", "gllpoints_all"),
                                            contiguous=False,
                                            chunksizes=([self.len_time, self.chunksize]))
        
        # Create disp_z variable
        disp_z = nc_out.groups['Snapshots'].createVariable("disp_z", np.float32,
                                            ("snapshots", "gllpoints_all"),
                                            contiguous=False,
                                            chunksizes=([self.len_time, self.chunksize]))
        
        # Create disp_p variable - not needed for PZ
        if source_type != 'PZ':
            disp_p = nc_out.groups['Snapshots'].createVariable("disp_p", np.float32,
                                                ("snapshots", "gllpoints_all"),
                                            contiguous=False,
                                            chunksizes=([self.len_time, self.chunksize]))
            
        # Populate the displacement variables
        self._load_displacement_data_to_netCDF(nc_out, source_type, processes)
        
        # Calculate the source time function
        STF, dSTF = tools.calculate_STF(
                self.source_info, self.element_group_info['metadata']['data_time'],
                self.source_magnitude)
        
        # Create the two STF variables
        stf_dump = nc_out.groups['Snapshots'].createVariable("stf_dump", np.float32,
                                            ("snapshots"))
        stf_dump[:] = STF
        stf_d_dump = nc_out.groups['Snapshots'].createVariable("stf_d_dump", np.float32,
                                            ("snapshots"))
        stf_d_dump[:] = dSTF



    def _load_displacement_data_to_netCDF(self, nc_out, source_type, processes):
        """
        This function loads the displacement data from the element output and
        assigns it to the netCDF file.
        
        Currently individual channels are loaded as a whole and assigned.
        Memory saving measures may be implemented in the future.
        
        This is done in parallel.
        """
        
        # Which Na is needed for each component.
        # These are the indices of the only displacement slice which is needed.
        if source_type == 'PZ':
            Nas = {'Us':0, 'Up':0, 'Uz':0}
        elif source_type == 'PX0':
            Nas = {'Us':1, 'Up':2, 'Uz':1}
        elif source_type == 'PX1':
            Nas = {'Us':2, 'Up':1, 'Uz':2}
        
        # Channels
        channels = self.element_group_info['metadata']['detailed_channels']
        
        # Na indices to access
        Na_indices = [Nas[channel] for channel in channels]

        # Take only unique GLL points
        coordinates, unique_indices, unique_inverse = self._get_unique_points(
                            return_all=True)

        # The files
        nc_files = self.element_group_info['metadata']['files']

        # Break files into lists for reading in parallel
        nc_file_lists = [nc_files[i::processes] for i in range(processes)]
        
        # Indices for each file to assign to
        element_length = len(self.element_group_info['inplane']['sem_inds'])
        index_limits = np.append(0, np.cumsum([f.sizes['dim_element'] * element_length
                                             for f in nc_files]))
        index_assignments = [[unique_inverse[np.arange(index_limits[i], 
                                                      index_limits[i+1])] 
                                    for i in range(len(index_limits)-1)][i::processes]
                                    for i in range(processes)]
        
        # Take only processes that have files allocated to them
        num_processes = np.where(np.array(
                            [len(x) for x in nc_file_lists]) != 0.)[0][-1] + 1
        nc_file_lists = nc_file_lists[:num_processes]
        index_assignments = index_assignments[:num_processes]
        
        # Loop through channels
        for ch_i, (channel, Na_i) in enumerate(zip(channels, Na_indices)):
            
            # No need for Up if source is PZ
            if not (channel == 'Up' and source_type == 'PZ'):
            
                # netCDF variable
                nc_key = "disp_" + channel[-1]
                var_out = nc_out.groups['Snapshots'][nc_key]
                
                # Load the displacement
                self._load_displacement_data_per_channel(nc_file_lists,
                                index_assignments, var_out, Na_i, ch_i, source_type)


        
    def _load_displacement_data_per_channel(self, nc_file_lists, index_assignments, 
                                                        var_out, Na_i, ch_i, source_type):
        """
        This function loads the displacement from a file list for a particular channel.
        """

        # Displacement array to assign to 
        displacement = np.empty((self.len_time, self.npoints))
        
        # Loop through files and read in parallel
        # Start a Queue for managing the outputs
        q = Queue()
        jobs = []
        for n, (nc_file_list, index_assignment_list) in enumerate(zip(
                                            nc_file_lists, index_assignments)):
            logger.info(f"   {var_out.name}: starting process number {n+1}, " + 
                        f"which has {len(nc_file_list)} files to read")
            p = Process(name=str(n+1), target=self._read_displacement_from_files,
                        args=(var_out.name, nc_file_list, ch_i, Na_i, 
                                    index_assignment_list, q))
            jobs.append(p)
            p.start()
            
        # For each process assign the output
        for n in range(len(jobs)):
            
            # Assign the output to the displacement array
            output, indices = q.get()
            displacement[:,indices] = output
            
        # Join each process to ensure we don't continue until all are finished
        for p in jobs:
            p.join()
        logger.info(f'   {var_out.name}: all processes have now finished')
        
        # Assign to netCDF output
        if source_type == 'PZ':
            var_out[:,:] = displacement
        elif 'PX' in source_type:
            var_out[:,:] = displacement * np.sqrt(np.pi)
        
        logger.info(f'  {var_out.name}: finished loading displacement')



    def _read_displacement_from_files(self, name, nc_file_list, ch_i, Na_i, 
                                            index_assignments_list, queue):
        """
        This function loads the displacement from the files in a file list.
        """
        
        # Create an output to assign to
        all_channel_displacements = [np.empty((x.shape[0], self.len_time),
                            dtype=np.float32) for x in index_assignments_list]
        
        # Loop through assigned files
        for i, nc_file in enumerate(nc_file_list):
        
            all_channel_displacements[i] = (
                np.array(nc_file['data_wave__NaG=%d' % self.Nr][
                                    :, Na_i, :, ch_i, :], dtype=np.float32)
                                ).reshape(-1, self.len_time)
        
        # Concatenate into one array
        channel_displacements = np.concatenate(all_channel_displacements, axis=0)
        
        # What indices to assign these to
        indices = np.concatenate(index_assignments_list)

        # Order by index, and transpose
        index_sorting = indices.argsort()
        channel_displacements = channel_displacements[index_sorting, :].T
        indices = np.sort(indices)
        
        # Add to the queue
        queue.put((channel_displacements, indices))
        logger.info(f"    {name}: process number " +
              f" {multiprocessing.current_process().name} has " + 
              f"finished having loaded data from {channel_displacements.shape[1]}" +
               " GLL points")



    def _create_output_mesh(self, nc_out):
        """
        Creates the `Mesh` group which contains, yes you guessed it, the information
        about the mesh.
        """
        
        # Load the mesh
        full_mesh = mesh.Exodus_Mesh(self.mesh_file,
                        self.element_group_info['metadata']['list_element'], 
                        self._get_element_mapping('Mids'),
                        self._get_element_mapping('FEM'),
                        self._get_element_mapping('SEM'),
                        self._get_unique_points(return_all=False),
                        basic_only=False,
                        output_everything=self.output_all_mesh_parameters)
        
        # Create the group
        nc_out.createGroup('Mesh')
        
        # Make the Mesh group dimensions
        # Number of elements in mesh
        nc_out.groups['Mesh'].createDimension('elements', self.nelements)
        # Control points per elements - for these meshes should always be the four nodes
        nc_out.groups['Mesh'].createDimension('control_points', 4)
        # Polynomial order plus 1
        nc_out.groups['Mesh'].createDimension('npol', 
                                        self.element_group_info['inplane']['npol'] + 1)
        
        # Calculate the GLL and GLJ points
        gll_points = tools.calculate_GLL(
                            self.element_group_info['inplane']['npol'] + 1)
        glj_points = tools.calculate_GLJ(
                            self.element_group_info['inplane']['npol'])
        
        # Next calculate the drivatives of the basis functions
        # Needed for the strain
        G2 = tools.calculate_G2(
                    self.element_group_info['inplane']['npol'], gll_points)
        G1 = tools.calculate_G1(
                    self.element_group_info['inplane']['npol'], glj_points)
        G0 = G1[0,:]
        
        # Make a dictionary of required properties and values
        properties = {"midpoint_mesh": {"dimensions": ("elements"),
                                        "dtype": np.int32,
                                        "key": "midpoint_mesh"},
                      "eltype": {"dimensions": ("elements"),
                                 "dtype": np.int32,
                                 "key": "element_type"},
                      "axis": {"dimensions": ("elements"),
                               "dtype": np.int32,
                                "key": "axis"},
                      "fem_mesh": {"dimensions": ("elements", "control_points"),
                                   "dtype": np.int32,
                                   "key": "fem_mesh"},
                      "sem_mesh": {"dimensions": ("elements", "npol", "npol"),
                                   "dtype": np.int32,
                                   "key": "sem_mesh"},
                      "mp_mesh_S": {"dimensions": ("elements"),
                                    "dtype": np.float32,
                                    "key": "mp_S"},
                      "mp_mesh_Z": {"dimensions": ("elements"),
                                    "dtype": np.float32,
                                    "key": "mp_Z"},
                      "G0": {"dimensions": ("npol"),
                             "dtype": np.float64,
                             "values": G0},
                      "G1": {"dimensions": ("npol", "npol"),
                             "dtype": np.float64,
                             "values": G1},
                      "G2": {"dimensions": ("npol", "npol"),
                             "dtype": np.float64,
                             "values": G2},
                      "glj": {"dimensions": ("npol"),
                              "dtype": np.float64,
                              "values": glj_points},
                      "gll": {"dimensions": ("npol"),
                              "dtype": np.float64,
                              "values": gll_points},
                      "mesh_S": {"dimensions": ("gllpoints_all"),
                                 "dtype": np.float32,
                                 "key": "S"},
                      "mesh_Z": {"dimensions": ("gllpoints_all"),
                                 "dtype": np.float32,
                                 "key": "Z"}}
        
        # Add seismic properties that are in the mesh
        # Note that not all were saved during mesh loading
        # Look in mesh.py for more details
        for param in full_mesh.values.keys():
            if param.isupper() and param not in ['S', 'Z']:
                properties['mesh_' + param.lower()] = {
                                 "dimensions": ("gllpoints_all"),
                                 "dtype": np.float32,
                                 "key": param}
        
        # Assign all properties to the netCDF file
        for prop in properties:
            variable = nc_out.groups['Mesh'].createVariable(
                                    prop , properties[prop]["dtype"],
                                    properties[prop]["dimensions"])
            if "values" in properties[prop]:
                variable[:] = properties[prop]["values"]
            elif "key" in properties[prop]:
                variable[:] = full_mesh.values[properties[prop]["key"]][:]



    def _create_database_attributes(self, source_type):
        """
        This function creates the root attributes for the Instaseis database.
        These are written by AxiSEM, but not all are required by Instaseis.
        They should nevertheless exist so as not to error.
        
        Attributes that I am unsure of are marked UNSURE with a guess at their meaning.
        """
        
        # Read in the mesh, but only the bare basics
        basic_mesh = mesh.Exodus_Mesh(self.mesh_file, basic_only=True)
        
        # Create an output to store things
        attributes = {}
        
        # Number of points in database
        attributes['npoints'] = self.npoints
        
        # Number of elements in mesh
        attributes['nelem_kwf_global'] = self.element_group_info[
                                    'metadata']['elements_index_limits'][-1]
        
        # DOES NOT MATTER
        attributes['file version'] = 10
        
        # Background model that was used
        attributes['background model'] = basic_mesh.name
        attributes['external model name'] = ''
        
        # Attenuation yes/no
        # 1 is true, 0 is false
        attributes['attenuation'] = int(basic_mesh.anelastic)
        
        # Radius of planet
        attributes['planet radius'] = basic_mesh.parameters['radius'] * 1e-3
        
        # Set the time - use now
        attributes['datetime'] = datetime.datetime.now(datetime.UTC).isoformat()
        
        # Boring computer things
        # These all DO NOT MATTER
        attributes['git commit hash'] = 'not important'
        attributes['user name'] = 'not important'
        attributes['host name'] = 'not important'
        attributes['compiler brand'] = 'not important'
        attributes['compiler version'] = 'not important'
        attributes['FFLAGS'] = 'not important'
        attributes['CFLAGS'] = 'not important'
        attributes['LDFLAGS'] = 'not important'
        attributes['OpenMP'] = 'not important'
        
        # Time scheme used for integration
        # Set to Newmark
        attributes['time scheme'] = 'newmark2'
        
        # Time steps
        dt_factor = float(self.element_group_info['temporal']['sampling_period'][3:])
        attributes['time step in sec'] = self.dt / dt_factor
        attributes['number of time steps'] = 0
        
        # Salvusmeshlite uses npol=4 as standard
        attributes['npol'] = self.element_group_info['inplane']['npol']
        
        # Excitation is always monopole and dipole for PZ and PX, respectively
        excitation_type = {'PZ': 'monopole',
                           'PX0': 'dipole',
                           'PX1': 'dipole'}
        attributes['excitation type'] = excitation_type[source_type]
        
        # Source type in AxiSEM language
        source_type_convert = {'PZ': 'vertforce',
                               'PX0': 'thetaforce',
                               'PX1': 'phiforce'}
        attributes['source type'] = source_type_convert[source_type]
        
        # Source time function
        # Convert to AxiSEM language
        source_conversion = {'ERF':'errorf',
                             'GAUSSIAN':'gauss_0',
                             'FIRST_DERIVATIVE':'gauss_1',
                             'RICKER':'gauss_2'}
        attributes['source time function'] = source_conversion[
                self.source_info['source_time_function']['use_derivative_integral']]
        
        # Simulation type should always be force
        attributes['simulation type'] = 'force'
        
        # Source period, use that of the mesh
        attributes['dominant source period'] = basic_mesh.parameters['minimum_period']
        
        # Source parameters - already checked so are fixed
        attributes['source depth in km'] = 0.
        attributes['Source colatitude'] = 0.
        attributes['Source longitude'] = 0.
        attributes['scalar source magnitude'] = self.source_magnitude
        
        # Don't bother delaring stations - not needed for Instaseis
        # DOES NOT MATTER
        attributes['number of receivers'] = 0
        attributes['length of seismogram  in time samples'] = attributes[
                                                            'number of time steps']
        attributes['seismogram sampling in sec'] = attributes['time step in sec']
        
        # Number of snapshots
        attributes['number of strain dumps'] = self.len_time
        
        # Time between snapshots
        attributes['strain dump sampling rate in sec'] = self.dt
        
        # Dump type is always displacement
        # IMPORTANT
        attributes['dump type (displ_only, displ_velo, fullfields)'] = 'displ_only'
        
        # Bounds of the wavefield - only one element group
        attributes['kernel wavefield rmin'] = self.element_group_info['elements'][
                                        'vertical_range'][0] * 1e-3
        attributes['kernel wavefield rmax'] = self.element_group_info['elements'][
                                        'vertical_range'][1] * 1e-3
        attributes['kernel wavefield colatmin'] = np.degrees(
                self.element_group_info['elements']['horizontal_range'][0])
        attributes['kernel wavefield colatmax'] = np.degrees(
                self.element_group_info['elements']['horizontal_range'][1])
        
        # UNSURE why this is zero
        attributes['number of snapshot dumps'] = 0
        attributes['snapshot dump sampling rate in sec'] = 0.
        
        # DOES NOT MATTER
        attributes['receiver components'] = 'cyl'
        
        # DOES NOT MATTER
        attributes['ibeg'] = 0
        attributes['iend'] = attributes['npol']
        attributes['jbeg'] = 0
        attributes['jend'] = attributes['npol']
        
        # Source shift parameters
        # This is essentially how long before the source time the database starts
        # IMPORTANT
        attributes['source shift factor in sec'] = (-1) * self.start_time
        attributes['source shift factor for deltat'] = (
                attributes['source shift factor in sec'] / 
                        attributes['time step in sec'])
        attributes['source shift factor for seis_dt'] = attributes[
                                                'source shift factor for deltat']
        attributes['source shift factor for deltat_coarse'] = (
                attributes['source shift factor in sec'] / 
                        attributes['strain dump sampling rate in sec'])
                        
        # Some more receiver things
        # DOES NOT MATTER
        attributes['receiver file type'] = 'stations'
        attributes['receiver spacing (0 if not even)'] = 0.
        
        # Some random things that must be true (or we wouldn't be doing this...)
        # DOES NOT MATTER
        attributes['use netcdf for wavefield output?'] = 'T'
        attributes['percent completed'] = 100
        attributes['finalized'] = 1
        
        return attributes
