#!/bin/bash
#SBATCH --job-name=MC
#SBATCH --account=pfaendtner
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --mem=150gb
# E-mail Notification, see man sbatch for options

## SBATCH --workdir=$SLURM_SUBMIT_DIR

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

module load ompi/4.1.3
mpif90 -ffree-form -fbounds-check -O3 f.f90 -o f.x

# the number of processes = the number of markov chains
mpirun -np 40 ./f.x


exit 0
