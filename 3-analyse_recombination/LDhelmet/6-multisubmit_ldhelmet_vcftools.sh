#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### Simple script for multiple submissions of vcftools to cluster #####

FILE_DIR="path/to/chr.vcf.gz/"

### In case you have forgotten to zip, this loop will do it for you.
#for filename in $FILE_DIR/chr*.vcf; do
#bgzip $filename ;done;

### Actual submission loop for vcftools --ldhelmet. See more in 6-ldhelmet_vcftools.sh
for zipped in $FILE_DIR/chr*.vcf.gz; do
sh ldhelmet_vcftools.sh $zipped ;done;
