#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00

##### Makes separate vcf files for each chromosome #####

directory="/path/to/directory/where/you/save/chr.vcfs/"

### For creating the files
grep -v "#" SwAsp_94samples.filter.gt.scaffolds.chr.filter.maf.hwe.vcf.recode.negative_strand_complemented.headered.vcf.25_random.recode.vcf | awk -v d=$directory 'BEGIN {FS="\t"};{close(f)};!seen[$1]++{f=d"/"$1".vcf"};{print $0 >> f}'

### For putting in the header
grep '#' SwAsp_94samples.filter.gt.scaffolds.chr.filter.maf.hwe.vcf.recode.negative_strand_complemented.headered.vcf.25_random.recode.vcf > header
for file in $directory/chr*.vcf; do
cat header $file > temp
mv temp $file ;done;
