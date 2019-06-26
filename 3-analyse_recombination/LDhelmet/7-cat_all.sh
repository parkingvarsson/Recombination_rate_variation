#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### This script loops through the chromosome folders you have the #####
##### haplotype fasta-files in and collects each file in these folders #####
##### into one fasta.all -file to be used for downstream #####

FILE_DIR="/path/to/haplotype/fastas/"

for filepath in $FILE_DIR/chr*; do
    cat $filepath/* > $filepath".all"
done
