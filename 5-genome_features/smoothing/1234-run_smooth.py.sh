#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### Script for running each of the smoothing scripts. They are intended #####
##### to be run separately. Comment the others out. #####

#python 1-smooth_recombination.py $1
#python 2-smooth_features.py $1
#python 3-smooth_methylation.py $1
python 4-smooth_substitutiondensity.py $1
