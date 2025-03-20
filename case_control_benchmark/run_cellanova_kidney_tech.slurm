#!/bin/bash

#SBATCH --partition=week

#SBATCH --job-name=case_cpu

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=6-00:00:00

#SBATCH --output=1_cellanova_kidney2.out



jupyter nbconvert --to notebook --execute cellanova_kidney2.ipynb
