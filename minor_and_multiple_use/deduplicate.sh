#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### Generic script for removing duplicates from a file #####

awk '!seen[$0]++' $1 > $1".deduplicated"
