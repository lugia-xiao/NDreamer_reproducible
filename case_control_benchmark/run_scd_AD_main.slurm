#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=scd_AD_main

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_scd_AD_main.out



jupyter nbconvert --to notebook --execute scd_main_AD.ipynb
