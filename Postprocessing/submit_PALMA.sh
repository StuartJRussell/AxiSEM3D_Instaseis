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
export time="01:00:00"

# Partition to submit to
# Recommended: bigsmp
export partition="express"

# Your email address
export email=""

##############################
#     --- Submission ---     #
# Do not change this section #
##############################

# Total number of tasks
totaltasks=$((nodnum * ntaskspernode))

# SLURM parameters
echo "#!/bin/bash -l"                                        > slurm_sub
echo "#SBATCH --job-name=AxiSEM3D_Conversion"               >> slurm_sub
echo "#SBATCH --nodes=$nodnum"                              >> slurm_sub
echo "#SBATCH --ntasks-per-node=$ntaskspernode"             >> slurm_sub
echo "#SBATCH --time=$time"                                 >> slurm_sub
echo "#SBATCH --no-requeue"                                 >> slurm_sub
echo "#SBATCH --partition=$partition"                       >> slurm_sub
echo "#SBATCH --mail-type=ALL"                              >> slurm_sub
echo "#SBATCH --mail-user=$email"                           >> slurm_sub

# Modules to load
echo "module --force purge"                                 >> slurm_sub
echo "module load palma/2024a"                              >> slurm_sub
echo "module load GCCcore/13.3.0"                           >> slurm_sub
echo "module load Python/3.12.3"                            >> slurm_sub

# This line
echo "export OMPI_MCA_pml=ucx"                              >> slurm_sub

# Run the thing
echo "python3 -u convert_to_Instaseis.py "$totaltasks" >& "$outname    >> slurm_sub

# Run the job
sbatch slurm_sub




