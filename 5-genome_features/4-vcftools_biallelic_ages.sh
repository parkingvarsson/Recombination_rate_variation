#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00

module load bioinfo-tools
module load vcftools/0.1.15

file_dir="/path/to/biallelic/vcf/"
tremuloides="tremuloides.list.txt"
for file in $file_dir/SRR1569781.HC-*; do
echo ${file%%.recode.vcf}.old
echo ${file%%.recode.vcf}.new
vcftools --vcf $file --positions $tremuloides --recode --out ${file%%.recode.vcf}.new
vcftools --vcf $file --exclude-positions $tremuloides --recode --out ${file%%.recode.vcf}.old
done
