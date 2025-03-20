#!/bin/bash



#SBATCH --job-name=AD_get_EMS

#SBATCH --output=1_AD_get_EMS.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=45G

#SBATCH --gpus=v100:1

#SBATCH --partition=gpu

#SBATCH --time=48:00:00



jupyter nbconvert --to notebook --execute AD_get_EMS2.ipynb
