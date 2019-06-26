#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 48:00:00

module load bioinfo-tools
module load vcftools/0.1.15

# $1 = vcf with polymorphic sites in P. tremula
# $2 = bed -file with neutral positions

vcftools --vcf $1 --SNPdensity 50000 --remove-indels --bed $2  --out ${1%.*}.neutral_intergenic
