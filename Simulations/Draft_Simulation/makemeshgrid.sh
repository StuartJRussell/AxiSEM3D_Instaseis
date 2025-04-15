#!/bin/bash

# Activate Conda environment
conda activate axisem3d_mesher

# Run mesher
python3 -m salvus_mesh_lite.interface --input_file=./input/inparam.mesh.yaml

# Move to input/
mv ./*.e input/

# Deactivate environment
conda deactivate
