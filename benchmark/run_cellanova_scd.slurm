#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=scd_cellanova_experimental

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_cellanova_scd.out



jupyter nbconvert --to notebook --execute scd.ipynb

jupyter nbconvert --to notebook --execute cellanova.ipynb
