#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=plot_batch_effect_case

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_plot_batch_effect.out



jupyter nbconvert --to notebook --execute plot_batch_effect.ipynb
