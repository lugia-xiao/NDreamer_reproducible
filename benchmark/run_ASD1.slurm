#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=ASD1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_ASD1.out



jupyter nbconvert --to notebook --execute run_ASD1.ipynb
