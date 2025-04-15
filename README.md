
# AxiSEM3D_Instaseis

This repository contains scripts that allow the element outputs of AxiSEM3D simulations to be converted to Instaseis database format. There are three directories:

- `Conda_Environments/` contains .yml files to create conda environments for this package, instaseis, and the mesher for AxiSEM3D.
- `Simulations/` contains the necessary scripts for setting up the simulations.
- `Postprocessing/` contains the scripts for converting the outputs of the simulations to Instaseis databases.

## Before using

Before using this package you must have a working AxiSEM3D executable and in order to use the simulation creation script this ust be located in `Simulations/Draft_Simulation`. The AxiSEM3D code and installation instructions can be found at https://github.com/AxiSEMunity/AxiSEM3D.

The AxiSEM3D executable and SLURM submission scripts that are included in the package are specific to the PALMA II cluster of the University of Muenster and should be replaced/edited. See the specific README files in the subdirectories for further guidance.

## Quick example

The following is a quick example of how to obtain an Instaseis database.

First ensure that the conda envirtonment is activated
```
conda activate axiinsta
```

Next navigate into the `Simulations/` directory and then edit the parameters in `make_simulation.py` and then run the script
```
python3 make_simulation.py
```

This should create a new directory which contains the AxiSEM3D simulations. Run this simulation, and once it is finished navigate into the `Postprocessing/` directory. Edit the parameters in `convert_to_Instaseis.py` and then run it. To run with 4 processes run 
```
python3 convert_to_Instaseis.py 4
```

This should create a `Databases/` directory in the root location, which contains the Instaseis databases.
