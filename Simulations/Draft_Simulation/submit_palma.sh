#!/bin/bash -l

###############################
#       --- Inputs ---        #
# You can change this section #
###############################

# Define an output file
export outname="OUTPUTFILE"

# Number of nodes and tasks per node
export nodnum=1
export ntaskspernode=36

# Time to allow
export time="00:30:00"

# Partition to submit to
export partition="express"

# Your email address
export email=""

##############################
#     --- Submission ---     #
# Do not change this section #
##############################

# Define an output file
export outputname="OUTPUTFILE"

######### write slurm submission file ##########

# Total number of tasks
totaltasks=$((nodnum * ntaskspernode))

# SLURM parameters
echo "#!/bin/bash -l"                                                > slurm_sub
echo "#SBATCH --job-name=AxiSEM3D"                                  >> slurm_sub
echo "#SBATCH --nodes=$nodnum"                                      >> slurm_sub
echo "#SBATCH --ntasks-per-node=$ntaskspernode"                     >> slurm_sub
echo "#SBATCH --time=$time"                                         >> slurm_sub
echo "#SBATCH --no-requeue"                                         >> slurm_sub
echo "#SBATCH --partition=$partition"                               >> slurm_sub
echo "#SBATCH --mail-type=ALL"                                      >> slurm_sub
echo "#SBATCH --mail-user=$email"                                   >> slurm_sub

# Modules to load
echo "module --force purge"                                         >> slurm_sub
echo "module load palma/2022a"                                      >> slurm_sub
echo "module load GCC/11.3.0"                                       >> slurm_sub
echo "module load CMake/3.23.1"                                     >> slurm_sub
echo "module load OpenMPI/4.1.4"                                    >> slurm_sub
echo "module load FFTW/3.3.10"                                      >> slurm_sub
echo "module load METIS/5.1.0"                                      >> slurm_sub
echo "module load netCDF/4.9.0"                                     >> slurm_sub

# This line
echo "export OMPI_MCA_pml=ucx"                                      >> slurm_sub

# Run the thing
echo "mpirun ./axisem3d_palma >& "$outputname                       >> slurm_sub

# Run the job
sbatch slurm_sub 
