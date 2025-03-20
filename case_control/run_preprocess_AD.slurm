#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=download_SEA-AD

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=16

#SBATCH --mem-per-cpu=128G

#SBATCH --time=24:00:00

#SBATCH --output=1_AD_preprocess.out



jupyter nbconvert --to notebook --execute SEA_AD_preprocess.ipynb
