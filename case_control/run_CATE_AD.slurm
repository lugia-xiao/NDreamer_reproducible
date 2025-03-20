#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=AD_estimate_CATE

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_AD_estimate_CATE.out



jupyter nbconvert --to notebook --execute AD_CATE.ipynb
