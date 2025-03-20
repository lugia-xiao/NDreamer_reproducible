#!/bin/bash



#SBATCH --job-name=reproducible

#SBATCH --output=1_gpu.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem-per-cpu=20G

#SBATCH --gpus=a5000:1

#SBATCH --partition=gpu

#SBATCH --time=48:00:00



#jupyter nbconvert --to notebook --execute PBMC_yale.ipynb
#jupyter nbconvert --to notebook --execute PBMC.ipynb
#jupyter nbconvert --to notebook --execute virus.ipynb
jupyter nbconvert --to notebook --execute ECCITE.ipynb
jupyter nbconvert --to notebook --execute ASD.ipynb
