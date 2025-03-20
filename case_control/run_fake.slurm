#!/bin/bash



#SBATCH --job-name=fake_case_control

#SBATCH --output=1_fake.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=16G

#SBATCH --gpus=v100:1

#SBATCH --partition=gpu

#SBATCH --time=48:00:00



jupyter nbconvert --to notebook --execute fake_signal_condition.ipynb

jupyter nbconvert --to notebook --execute fake_signal_all.ipynb
