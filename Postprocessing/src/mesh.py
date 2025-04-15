
# Import modules
import logging
import numpy as np
import netCDF4 as nc
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

# Logger and configuration
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class Exodus_Mesh:
    """
    Class for handling the mesh files for AxiSEM3D simulations when creating Instaseis
    compatible reciprocal Green's function databases.
    
    Assumes that your model is an Exodus mesh created by SalvusMeshLite and is stored
    in the /input/ directory.
    
    This is quite slow for large meshes, so I might parallelise it in the future.
    """

    def __init__(self, mesh_file, element_indices=None, midpoint_mapping=None,
                fem_mapping=None, sem_mapping=None, points=None, basic_only=False,
                output_everything=False):
        
        # Load in the mesh file
        self.exodus_file = mesh_file
        self.output_everything = output_everything
        
        # If basic only then do that
        if basic_only:
            self._read_exodus_mesh_basic()
        
        # Otherwise read full mesh
        else:
            
            # Error if the inputs are None
            checks = [element_indices is None,
                      midpoint_mapping is None,
                      fem_mapping is None,
                      sem_mapping is None,
                      points is None]
            if any(checks):
                raise ValueError("Some of the element parameters are None")
            
            # Read the mesh
            self._read_exodus_mesh_full(element_indices, midpoint_mapping, fem_mapping,
                                                sem_mapping, points)



    def _read_exodus_mesh_full(self, element_indices, midpoint_mapping, fem_mapping,
                                                sem_mapping, points):
        """
        Reads in exodus mesh file used for AxiSEM3D simulation.
        
        The mesh is only saved on the points asked for. Currently the values at the
        points are calculated by simple linear interpolation as a function of radius and
        theta..
        """
        
        logger.info(f'  Reading the mesh file')
        
        # Open the .e file
        mesh = nc.Dataset(self.exodus_file)
        
        # Some mesh variables are just not needed by Instaseis and can be omitted
        # in order to save memory - list these here. 
        # To include a variable just remove it from the list or to output all set
        # `output_everything' to True
        if not self.output_everything:
            not_needed = ['dt', 'edge_aspect_ratio', 'equiangular_skewness',
                          'fluid', 'minimum_period', 'vp_dt', 'QMU', 'QKAPPA',
                          'RHO', 'VP', 'VS', 'VPH', 'VPV', 'VSH', 'VSV',
                          'XI', 'PHI', 'ETA', 'LAMBDA']
        else:
            not_needed = []
        
        # Mesh variable names
        var_names = [''.join(x) for x in 
                    mesh['name_elem_var'][:].data.astype(str)]
        
        # Global variables
        glo_var_names = [''.join(x) for x in 
                    mesh['name_glo_var'][:].data.astype(str)]
        
        # Sort the mesh so that the element midpoint indices are monotonically 
        # increasing - this reduces Instaseis access time by ~1/3.
        sorting = midpoint_mapping.argsort()
        element_indices = element_indices[sorting]
        fem_mapping = fem_mapping[sorting]
        sem_mapping = sem_mapping[sorting]
        midpoint_mapping = np.sort(midpoint_mapping)
        
        # Coordinates of nodes
        # Assume that z is redundant
        mesh_locations_nodes = np.vstack(
                (mesh['coordx'][:].data, mesh['coordy'][:].data)
                    ).T
        
        # Coordinates of nodes per element
        mesh_locations_elements_nodes = mesh_locations_nodes[
                                    mesh['connect1'][:].data-1][element_indices]
        
        # Is each element along the axis
        # These are treated differently by Instaseis
        axis = np.array(np.any(mesh_locations_elements_nodes[:,:,0] < 1., axis=1),
                                                dtype=np.int64)
        
        # Midpoint of each element
        element_midpoints = points[midpoint_mapping]
        
        # Reshape to be all points serialy
        mesh_points = mesh_locations_elements_nodes.reshape(-1, 2)
        
        # Take as R and theta
        mesh_points_rth = np.vstack((np.sqrt(np.sum(mesh_points ** 2, axis=1)),
                  np.degrees(np.arctan2(mesh_points[:,0], mesh_points[:,1])))).T
        
        # Points as R and theta
        points_rth = np.vstack((np.sqrt(np.sum(points ** 2, axis=1)),
                  np.degrees(np.arctan2(points[:,0], points[:,1])))).T
        
        # Seismic parameters are defined in four separate arrays, one for each node
        # so create an array for each element
        params = list(set([x.split('_')[0] for x in var_names if x.isupper()]))
        var_keys = [[f'vals_elem_var{i+1}eb1' for i, x in enumerate(var_names) 
                            if par in x] for par in params]
        
        # Get the values of the parameters on the mesh
        # Only take ones that are required by Instaseis
        logger.info(f'    Interpolating mesh parameters')
        values = {}
        for param, var_key in zip(params, var_keys):
            values_in = np.array([mesh[key][:].data[0, element_indices] 
                                            for key in var_key]).T.reshape(-1)
            values[param] = self._interpolate_parameters(
                                           points_in=mesh_points_rth,
                                           values_in=values_in,
                                           points_out=points_rth)
        
        # Set whether the mesh is anisotropic or not based on parameters
        self.anisotropic = 'ETA' in values
        
        # Set whether the mesh is anelastic or not based on parameters
        self.anelastic = 'QMU' in values
        
        # Calculate different properties only if needed
        # AxiSEM seems to just use VPH and VSH for me and lambda in get_model.F90
        if 'MU' not in not_needed:
            vs_var = 'VSH' if self.anisotropic else 'VS'
            vp_var = 'VPH' if self.anisotropic else 'VP'
            values['MU'] = values['RHO'][:] * values[vs_var][:] ** 2.
        if 'LAMBDA' not in not_needed:
            vs_var = 'VSH' if self.anisotropic else 'VS'
            vp_var = 'VPH' if self.anisotropic else 'VP'
            values['LAMBDA'] = (values['RHO'][:] * values[vp_var][:] ** 2.
                             - (values['RHO'][:] * values[vs_var][:] ** 2.))
        if 'XI' not in not_needed:
            if not self.anisotropic:
                values['XI'] = np.ones(values['VP'][:].shape)
            else:
                values['XI'] = values['VSH'][:] ** 2. / values['VSV'][:] ** 2.
        if 'PHI' not in not_needed:
            if not self.anisotropic:
                values['PHI'] = np.ones(values['VP'][:].shape)
            else:
                values['PHI'] = values['VPV'][:] ** 2. / values['VPH'][:] ** 2.
        if 'ETA' not in not_needed:
            if not self.anisotropic:
                values['ETA'] = np.ones(values['VP'][:].shape)
        
        # Assign whatever variables are needed for the mesh itself
        self.values = {}
        for param in values:
            if param not in not_needed:
                self.values[param] = values[param]
        
        # Assign the other variables
        other_vars = [x for x in var_names if x.split('_')[0] not in params]
        var_keys = [f'vals_elem_var{i+1}eb1' for i, x in enumerate(var_names) 
                            if x in other_vars]
        for param, var_key in zip(other_vars, var_keys):
            if param not in not_needed:
                self.values[param] = mesh[var_key][:].data[0, element_indices]
        
        # Assign coordinate variables
        self.values['S'] = points[:,0]
        self.values['Z'] = points[:,1]
        self.values['mp_S'] = element_midpoints[:,0]
        self.values['mp_Z'] = element_midpoints[:,1]
        self.values['axis'] = axis
        self.values['midpoint_mesh'] = midpoint_mapping
        self.values['fem_mesh'] = fem_mapping
        self.values['sem_mesh'] = sem_mapping
        
        # Read parameters
        self.parameters = {}
        for i, var in enumerate(glo_var_names):
            self.parameters[var] = mesh['vals_glo_var'][0, :].data[i]
        
        # Save mesh name
        self.name = getattr(mesh, 'title')
        
        logger.info(f'  Mesh file read')



    def _interpolate_parameters(self, points_in, values_in, points_out):
        """
        This function just gets interpolated values on the mesh at given points.
        
        This is in two stages as nans are sometimes returned at the edge of the mesh.
        """
        
        # Initialise interpolators
        interpolator_linear = LinearNDInterpolator(points_in, values_in)
        interpolator_nearest = NearestNDInterpolator(points_in, values_in)
        
        # Get the values at these points from linear interpolator
        values_out = interpolator_linear(points_out)
        
        # Some may be nan and need replacing, which are these
        nan_indices = np.where(np.isnan(values_out))
        
        # Replace them with nearest values
        values_out[nan_indices] = interpolator_nearest(points_out[nan_indices])
        
        return values_out



    def _read_exodus_mesh_basic(self):
        """
        Reads in exodus mesh file used for AxiSEM3D simulation.
        
        Only saves a small selection of essential parameters.
        """
        
        # Open the .e file
        mesh = nc.Dataset(self.exodus_file)
        
        # Save mesh name
        self.name = getattr(mesh, 'title')
        
        # Mesh variable names
        var_names = [''.join(x) for x in 
                    mesh['name_elem_var'][:].data.astype(str)]
        
        # Global variables
        glo_var_names = [''.join(x) for x in 
                    mesh['name_glo_var'][:].data.astype(str)]
        
        # Read parameters
        self.parameters = {}
        for i, var in enumerate(glo_var_names):
            self.parameters[var] = mesh['vals_glo_var'][0, :].data[i]
        
        # Set whether the mesh is anisotropic or not based on parameters
        self.anisotropic = 'ETA_0' in var_names
        
        # Set whether the mesh is anelastic or not based on parameters
        self.anelastic = 'QMU_0' in var_names
