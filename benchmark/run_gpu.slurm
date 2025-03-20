#!/bin/bash

#SBATCH --job-name=benchmark_gpu

#SBATCH --output=1_benchmark_gpu.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --gpus=a5000:1

#SBATCH --partition=gpu

#SBATCH --time=24:00:00



jupyter nbconvert --to notebook --execute scCAPE.ipynb

jupyter nbconvert --to notebook --execute scGen.ipynb
