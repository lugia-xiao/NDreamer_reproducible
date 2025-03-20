#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=evaluate_cellanova_experiment

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_evaluate_cellanova.out


jupyter nbconvert --to notebook --execute evaluate_cellanova.ipynb
