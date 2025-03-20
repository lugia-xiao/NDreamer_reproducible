#!/bin/bash

#SBATCH --partition=week

#SBATCH --job-name=scd_get_space_case_control

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=3-00:00:00

#SBATCH --output=1_scd_get_space.out



jupyter nbconvert --to notebook --execute scd_get_latent_case_control.ipynb
