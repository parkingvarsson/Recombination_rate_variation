#!/bin/bash -l
#SBATCH -A 
#SBATCH --mem=32G
#SBATCH -c 32

module load R

Rscript run-1.R $1


