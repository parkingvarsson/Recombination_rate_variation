#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### General script for making a bed out of any name and position #####
##### containing file #####

awk '{print $1"\t"$2"\t"$2+1}' $1 > ${1%%.txt}.bed
