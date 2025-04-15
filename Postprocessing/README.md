
# Postprocessing Directory

It is in this directory that all the scripts to convert the AxiSEM3D output to instaseis compatible databases are contained. It is as simple as editing the simulation name and path in `convert_to_Instaseis.py` and then running this script.

Locally this can be run by:
```
python3 convert_to_Instaseis.py 4
```
where the final number is the number of parallel processes to run.

Alternatively for very large simulations it will be necessary to also run this step on a HPC; `submit_PALMA.sh` is an example SLURM submission script that is compatible with the PALMA II cluster of the University of Muenster. It will need to be edited/replaced for use with other clusters.

The output databases will be contained in the `../Databases/` directory.
