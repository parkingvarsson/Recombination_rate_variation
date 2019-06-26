#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### Script for producing the vcf with 25 randomly chosen individuals. #####
##### You will find our randomly chosen individuals in the input_files: #####
##### SwAsp_94samples.filter.gt.scaffolds.chr.filter.maf.hwe.vcf.recode.negative_strand_complemented.headered.vcf.25_random.recode.vcf #####

module load bioinfo-tools
module load vcftools/0.1.15

vcftools --vcf "/path/to/complemented/vcf/file" --max-indv 25 --recode --out anyname".25_random"
