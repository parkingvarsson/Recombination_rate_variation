#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### Simple filtering for vcf without indels or one with, but indels need #####
##### to be removed sooner or later #####

vcftools --vcf $1 --maf 0.05 --hwe 0.002  --recode --out $1".maf.hwe"
