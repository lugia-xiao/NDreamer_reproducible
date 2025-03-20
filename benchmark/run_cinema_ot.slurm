#!/bin/bash



#SBATCH --partition=bigmem



#SBATCH --job-name=cinema_ot



#SBATCH --ntasks=1



#SBATCH --cpus-per-task=4



#SBATCH --mem-per-cpu=64G



#SBATCH --time=24:00:00



#SBATCH --output=1_cinema_ot_ITE.out







jupyter nbconvert --to notebook --execute cinema_ot_ITE.ipynb
