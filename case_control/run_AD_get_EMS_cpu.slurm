#!/bin/bash

#SBATCH --partition=week

#SBATCH --job-name=AD_decompse

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=6-00:00:00

#SBATCH --output=1_AD_decompse.out



jupyter nbconvert --to notebook --execute AD_estimate_ITE.ipynb
