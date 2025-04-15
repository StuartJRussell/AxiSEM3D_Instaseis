
# AxiSEM3D_Instaseis

This repository contains scripts that allow the element outputs of AxiSEM3D simulations to be converted to Instaseis database format. There are three directories:

- `Conda_Environments/` contains .yml files to create conda environments for this package, instaseis, and the mesher for AxiSEM3D.
- `Simulations/` contains the necessary scripts for setting up the simulations.
- `Postprocessing/` contains the scripts for converting the outputs of the simulations to Instaseis databases.

## Quick example

The following is a quick example of how to obtain an Instaseis database.

First ensure that the conda envirtonment is activated
```
conda activate axiinsta
```

Next navigate into the simulations directory
```
cd Simulations
```

Then edit the parameters in `make_simulation.py` and then run the script
```
python3 make_simulation.py
```

This should create a new directory which contains the AxiSEM3D simulations. Run this simulation.
