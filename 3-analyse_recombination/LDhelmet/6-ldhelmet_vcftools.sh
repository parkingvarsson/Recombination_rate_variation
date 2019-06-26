#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00

module load bioinfo-tools
module load vcftools/0.1.15

##### Script for producing snps and pos -files with vcftools #####

vcf=$1 # $1 is the gzipped vcf file used

### Extracts information from the filename
vcfname="${vcf##*/}"
chrom="${vcfname%%.*}"
chr="${chrom%?}"

### Performs vcftools --ldhelmet
vcftools --gzvcf $vcf --chr $chr --ldhelmet --out $chrom
