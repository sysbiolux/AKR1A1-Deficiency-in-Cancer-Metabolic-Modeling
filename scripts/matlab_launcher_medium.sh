#!/bin/bash -l
#SBATCH -J sampling_medium
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=evelyn.gonzalez@uni.lu
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --time=48:00:00
#SBATCH -p batch

# Load the module MATLAB
module purge
resif-load-swset-legacy  # load older HPC software versions

module load math/MATLAB
module load math/CPLEX #/12.10-GCCcore-9.3.0-Python-3.8.2 # loading cplex

export ILOG_CPLEX_PATH='/opt/apps/resif/aion/2020a/epyc/software/CPLEX'

# second form with CLI options '-r <input>' and '-logfile <output>.out', this can be used if 'quit' is used at the end of the matlab script
#srun -c ${SLURM_CPUS_PER_TASK} matlab -nodisplay -r ${HOME}/matlab_parallel/sampling_test_medium -logfile ${SCRATCH}/output_job_${SLURM_JOB_ID}.out

# first run with operators '< script.m > <output>.out'
srun -c ${SLURM_CPUS_PER_TASK} matlab -nodisplay -nosplash < ${HOME}/sampling_test_jeff/sampling_test_medium.m > ${SCRATCH}/output_job_${SLURM_JOB_ID}.out

# safeguard (!) afterwards
rm -rf ${HOME}/.matlab
rm -rf ${HOME}/java*
