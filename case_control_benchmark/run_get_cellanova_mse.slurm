#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=get_cellanova_mse_case_control

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_cellanova_mse.out



jupyter nbconvert --to notebook --execute get_cell_anova_mse.ipynb
