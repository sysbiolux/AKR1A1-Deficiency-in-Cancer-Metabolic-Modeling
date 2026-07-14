#!/bin/bash -l
#SBATCH -J sampling_medium
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=evelyn.gonzalez@uni.lu
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --time=48:00:00
#SBATCH -p batch

# Load software
module purge
resif-load-swset-legacy

module load math/MATLAB
module load math/CPLEX

export ILOG_CPLEX_PATH=/opt/apps/resif/aion/2020a/epyc/software/CPLEX

# Run MATLAB sampling script (standard models)
srun -c ${SLURM_CPUS_PER_TASK} \
matlab -nodisplay -nosplash \
-r "run('scripts/sampling_medium.m');"

# Cleanup
rm -rf ${HOME}/.matlab
rm -rf ${HOME}/java*