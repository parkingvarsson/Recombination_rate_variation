#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00

##### Simple script for making bed files from vcfs #####

module load bioinfo-tools
module load BEDOPS/2.4.3

vcf2bed --do-not-split < $1
