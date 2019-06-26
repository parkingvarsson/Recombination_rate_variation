#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00

module load bioinfo-tools
module load vcftools/0.1.15

vcftools --vcf $1 --SNPdensity 50000 --out ${1%%.recode.vcf}
