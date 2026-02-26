#!/bin/bash 
#SBATCH --job-name=gpu-smoke 
#SBATCH --output=gpu-smoke-%j.out 
#SBATCH --error=gpu-smoke-%j.err 
#SBATCH --partition=gpu 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --gpus-per-node=1 
#SBATCH --time=00:05:00


hostname echo 
"===== nvidia-smi -L =====" 
nvidia-smi -L 


source /etc/profile.d/west.sh
