#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00

for file in /path/to/cov/files/
    do
        sbatch 2-liftOver.sh $file
    done
