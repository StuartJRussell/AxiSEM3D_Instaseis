
# Conda Environments Directory

This directory contains files to setup conda environments that are required for this package to work.

The conda environments can be created from the .yml files by running the folling line, substituting the correct filename.

## axi_env.yml
Conda environment for running the scripts of this package. Install with:
```
conda env create -f axi_env.yml
```

## instaseis_env.yml
Conda environment that works with the instaseis Python module. This is useful for accessing the created databases. Install with:
```
conda env create -f instaseis_env.yml
```

## mesher_env.yml
Conda environment for the mesher for AxiSEM3D. This is required to created AxiSEM3D simulations. Install with:
```
conda env create -f mesher_env.yml
```

It is also necessary to install the mesher itself. This can be done by running the following lines:
```
conda activate axisem3d_mesher
pip install https://gitlab.com/Salvus/SalvusMeshLite/-/archive/master/SalvusMeshLite-master.zip
conda deactivate
```
