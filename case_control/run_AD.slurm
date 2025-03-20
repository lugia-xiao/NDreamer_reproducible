#!/bin/bash



#SBATCH --job-name=AD_case_control_a5000

#SBATCH --output=3.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=24G

#SBATCH --gpus=a5000:1

#SBATCH --partition=gpu

#SBATCH --time=2-00:00:00



jupyter nbconvert --to notebook --execute AD.ipynb
