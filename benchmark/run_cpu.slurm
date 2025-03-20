#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=benchmark_cpu

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1_benchmark_cpu.out



jupyter nbconvert --to notebook --execute cinema_ot.ipynb

jupyter nbconvert --to notebook --execute mixscape.ipynb

jupyter nbconvert --to notebook --execute scCAPE.ipynb

jupyter nbconvert --to notebook --execute scGen.ipynb
