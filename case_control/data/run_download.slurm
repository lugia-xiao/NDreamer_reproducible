#!/bin/bash

#SBATCH --partition=bigmem

#SBATCH --job-name=download_SEA-AD

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=64G

#SBATCH --time=24:00:00

#SBATCH --output=1.out



wget https://datasets.cellxgene.cziscience.com/c32964d2-3339-441f-8e56-7177234c7876.h5ad
