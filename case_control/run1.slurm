#!/bin/bash



#SBATCH --job-name=1case_control

#SBATCH --output=2.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=20G

#SBATCH --gpus=a5000:1

#SBATCH --partition=gpu

#SBATCH --time=48:00:00



jupyter nbconvert --to notebook --execute Mouse.ipynb
