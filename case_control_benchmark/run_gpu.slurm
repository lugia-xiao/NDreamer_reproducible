#!/bin/bash



#SBATCH --job-name=case_gpu

#SBATCH --output=2.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=32

#SBATCH --mem-per-cpu=25G

#SBATCH --gpus=a5000:1

#SBATCH --partition=scavenge_gpu

#SBATCH --time=24:00:00



jupyter nbconvert --to notebook --execute scdisinfact.ipynb
