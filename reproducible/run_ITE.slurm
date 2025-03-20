#!/bin/bash

#SBATCH --partition=week

#SBATCH --job-name=evaluate_ITE

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=6-00:00:00

#SBATCH --output=1_ITE.out



jupyter nbconvert --to notebook --execute evaluate_ITE.ipynb
