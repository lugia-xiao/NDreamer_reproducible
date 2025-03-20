#!/bin/bash



#SBATCH --job-name=AD_case_control_A100

#SBATCH --output=4.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=24G

#SBATCH --gpus=a100:1

#SBATCH --partition=gpu

#SBATCH --time=2-00:00:00



jupyter nbconvert --to notebook --execute AD_ADNC.ipynb
