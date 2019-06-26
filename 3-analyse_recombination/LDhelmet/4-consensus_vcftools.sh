#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00

module load bioinfo-tools
module load vcftools/0.1.15

##### This script runs the vcftools vcf-consensus -script to gain #####
##### haplotype fasta files for each of the 19 chromosomes in each of #####
##### 25 random individuals. This means there will be 50 files generated #####
##### for each chromosome. The looper for this exists in
##### 4-bgzip_and_tabix.sh #####

# $1 should be the reference fasta-file
# $2 should be the vcf file with the 25 random individuals
# $3 is the name of the 1/25 random individuals currently processed

### Extracts the chromosome name from the vcf and creates folder for it
vcfname="${2##*/}"
chrom="${vcfname%%.*}"

for chromosome in $chrom; do
mkdir $chrom ;done;

>$chrom.$sample.25_rand.fa

### Performs vcf-consensus for both haplotypes and adjusts the names in the
### output files
cat $1 | vcf-consensus --haplotype 1 --sample $3 $2 >> $chrom.$sample.25_rand.fa
sed -i -e "s/$chrom/$sample.1/g" $chrom.$sample.25_rand.fa
cat $1 | vcf-consensus --haplotype 2 --sample $3 $2 >> $chrom.$sample.25_rand.fa
sed -i -e "s/$chrom/$sample.2/g" $chrom.$sample.25_rand.fa

### Note ###
# --haplotype 1==Reference
# --haplotype 2==Alternate
