#!/bin/bash

#SBATCH --partition=week

#SBATCH --job-name=evaluate_case_benchmark

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=6

#SBATCH --mem-per-cpu=64G

#SBATCH --time=6-00:00:00

#SBATCH --output=1_evaluate.out



jupyter nbconvert --to notebook --execute evaluate_cellanova.ipynb
