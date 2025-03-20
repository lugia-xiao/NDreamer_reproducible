#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=plot_benchmark

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_plot_AD.out



jupyter nbconvert --to notebook --execute plot_cellanova_AD.ipynb

jupyter nbconvert --to notebook --execute plot_cellanova.ipynb

jupyter nbconvert --to notebook --execute plot_scd_AD.ipynb

jupyter nbconvert --to notebook --execute plot_scd.ipynb
