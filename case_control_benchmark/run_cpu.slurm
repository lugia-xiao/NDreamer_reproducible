#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=case_cpu

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=14

#SBATCH --mem-per-cpu=140G

#SBATCH --time=24:00:00

#SBATCH --output=1.out



jupyter nbconvert --to notebook --execute cellanova.ipynb
