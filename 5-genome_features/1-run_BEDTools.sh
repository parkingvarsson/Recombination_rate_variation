#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00

module load bioinfo-tools
module load BEDTools/2.27.1

##### This script is used to produce gene and repeat coverages for windows #####

bedtools makewindows -g $1 -w 50000  >  windows.bed

dir="/path/to/gffs/"

for gff in $dir/*gff3; do
    bedtools annotate -i windows.bed -files $gff > $gff.coverage.bed
    #echo $gff
done
