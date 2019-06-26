#!/bin/bash -l

#SBATCH -c 1
#SBATCH -A 
#SBATCH --mail-user 
#SBATCH --mail-type=ALL

set -x

module load R

LG=$1

Rscript Consensus_Map.R ${LG}
