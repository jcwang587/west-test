#!/bin/bash 

#SBATCH --job-name=gpu-smoke 
#SBATCH --output=gpu-smoke-%j.out 
#SBATCH --error=gpu-smoke-%j.err 
#SBATCH --partition=gpu 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --gpus-per-node=1 
#SBATCH --time=00:05:00

########################################################
# WEST SIMULATION SETUP
########################################################

# Print GPU node information
nvidia-smi -L 

# Initialize WEST environment
source /etc/profile.d/west.sh

# Set CUDA visible devices and OpenMP threads
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1

# Run WEST simulation
mpirun -np 1 west.x -i west.in > west.out 2> west.err



