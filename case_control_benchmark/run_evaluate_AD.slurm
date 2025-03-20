#!/bin/bash

#SBATCH --partition=week

#SBATCH --job-name=evaluate_case_benchmark_AD

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=6

#SBATCH --mem-per-cpu=80G

#SBATCH --time=6-00:00:00

#SBATCH --output=1_evaluate_AD.out



jupyter nbconvert --to notebook --execute evaluate_cellanova_AD.ipynb
