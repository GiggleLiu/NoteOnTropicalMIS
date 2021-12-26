#!/bin/sh

#SBATCH -p gpu
#SBATCH --gres=gpu:4
#SBATCH -c 5
#SBATCH --mem 16000
#SBATCH -t 0-1:00 # time (D-HH:MM)
#SBATCH --mail-user=jinguoliu@g.harvard.edu
module load cuda/11.4.2-fasrc01
julia -e "println(1)"
JULIA_NUM_THREADS=5 JULIA_CUDA_USE_BINARYBUILDER=false julia --project=.. sequencing_multigpu.jl true

