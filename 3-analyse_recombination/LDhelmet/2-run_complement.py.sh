#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### See more information in "2-complement.py" #####

python 2-complement.py $1 $2 # $1 = vcf to complement $2 = agp to compement by
