#!/bin/bash



#SBATCH --job-name=reproducible_ASD1



#SBATCH --output=1_ASD1.out



#SBATCH --ntasks=1



#SBATCH --cpus-per-task=4



#SBATCH --mem-per-cpu=20G



#SBATCH --gpus=v100:1



#SBATCH --partition=gpu



#SBATCH --time=48:00:00





jupyter nbconvert --to notebook --execute ASD1_nobatch.ipynb
