#!/bin/bash

#SBATCH --partition=day

#SBATCH --job-name=evaluate_exp_batch

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_evaluate_batch.out


jupyter nbconvert --to notebook --execute evaluate_batch.ipynb
