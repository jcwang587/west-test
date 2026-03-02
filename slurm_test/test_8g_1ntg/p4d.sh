#!/bin/bash 

#SBATCH --job-name=gpu-smoke 
#SBATCH --output=gpu-smoke-%j.out 
#SBATCH --error=gpu-smoke-%j.err 
#SBATCH --partition=gpu 
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

########################################################
# WEST SIMULATION SETUP
########################################################

# Print GPU node information
nvidia-smi -L 

# Initialize WEST environment
source /etc/profile.d/west.sh

# Set CUDA visible devices and OpenMP threads
export OMP_NUM_THREADS=1

# Run WEST simulation
mpirun -np 8 pw.x -pd true -ntg 1 -i pw.in > pw.out 2> pw.err